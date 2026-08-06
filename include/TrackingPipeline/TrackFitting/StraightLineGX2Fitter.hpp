#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include <Acts/Geometry/GeometryIdentifier.hpp>

#include <cstddef>
#include <cstdint>
#include <unordered_map>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/detail/BinningValueUtils.hpp"
#include "TrackingPipeline/Utilities/TupleHash.hpp"

/// @brief simple global chi2 fitter for intermediate estimations
///
/// The class performs a simple chi2 fit of a particle track accounting for
/// multiple scattering effects. The fitter relies on the initial estimate of
/// the track's direction.
class StraightLineGX2Fitter {
 public:
  /// @brief config struct
  struct Config {
    /// Coordinate indices
    std::size_t primaryIdx;
    std::size_t longIdx;
    std::size_t shortIdx;
    /// Geometry ID of the fit starting surface
    Acts::GeometryIdentifier startSurfaceGeoId;
    /// Surface map
    std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
        surfaceMap;
  };

  /// Enum with the GX2 estimate dimension indices
  enum GX2EstimateIndices : std::size_t {
    ePos0 = 0,
    ePos1 = 1,
    eDir0 = 2,
    eDir1 = 3,
    eEstimateSize = 4
  };

  /// Enum with the free dimension indices
  enum GX2FreeIndices : std::size_t {
    eFreePos0 = 0,
    eFreePos1 = 1,
    eFreePos2 = 2,
    eFreeDir0 = 3,
    eFreeDir1 = 4,
    eFreeDir2 = 5,
    eFreeSize = 6
  };

  /// Fitter covariance shorthand
  using Covariance = Acts::ActsSquareMatrix<eFreeSize>;

  /// @brief constructor
  ///
  /// @param cfg configuration struct
  explicit StraightLineGX2Fitter(const Config& cfg) : m_cfg(cfg) {
    for (const auto& [geoId, surf] : m_cfg.surfaceMap) {
      m_surfaces.emplace_back(geoId, surf);
    }
    std::sort(m_surfaces.begin(), m_surfaces.end(),
              [](const auto& a, const auto& b) {
                return (a.first.sensitive() < b.first.sensitive());
              });
  }

  /// @brief perform global chi2 fit
  ///
  /// @tparam container_t source link container type
  ///
  /// @param gctx geometry context
  /// @param sourceLinksContainer source links container
  /// @param thetaMcsRms rms scattering angle in the layers
  /// @param pos track position estimate to be updated
  /// @param dir track direction estimate to be updated
  /// @param cov track covariance estimate to be updated
  ///
  /// @return chi2 of the track fit
  template <typename container_t>
  double gx2Fit(const Acts::GeometryContext& gctx,
                const container_t& sourceLinkContainer, double thetaMcsRms,
                Acts::Vector3& pos, Acts::Vector3& dir, Covariance& cov) {
    Acts::Vector3 primaryDir = detail::indexToDirection(m_cfg.primaryIdx);

    std::unordered_map<Acts::GeometryIdentifier, Acts::Vector3> surfaceCenters;
    std::unordered_map<Acts::GeometryIdentifier, Acts::Vector3> surfaceNormals;
    for (const auto& [geoId, surf] : m_cfg.surfaceMap) {
      surfaceCenters.insert({geoId, surf->center(gctx)});
      surfaceNormals.insert(
          {geoId, surf->normal(gctx, surfaceCenters.at(geoId), primaryDir)});
    }

    std::unordered_map<std::pair<std::uint32_t, std::uint32_t>, double,
                       detail::PairHash>
        interSurfaceDistanceSums;
    std::size_t nLayers = m_surfaces.size();
    for (std::size_t i = 0; i < nLayers; i++) {
      const auto& [geoIdI, surfI] = m_surfaces.at(i);
      for (std::size_t j = 0; j < nLayers; j++) {
        const auto& [geoIdJ, surfJ] = m_surfaces.at(j);
        double sum = 0;
        for (std::size_t k = 0; k < std::min(i, j); k++) {
          const auto& [geoIdK, surfK] = m_surfaces.at(k);
          double dki = (surfaceCenters.at(geoIdI) - surfaceCenters.at(geoIdK))
                           .dot(surfaceNormals.at(geoIdI)) /
                       dir.dot(surfaceNormals.at(geoIdI));
          double dkj = (surfaceCenters.at(geoIdJ) - surfaceCenters.at(geoIdK))
                           .dot(surfaceNormals.at(geoIdJ)) /
                       dir.dot(surfaceNormals.at(geoIdJ));
          sum += dki * dkj;
        }
        interSurfaceDistanceSums.insert(
            {{geoIdI.sensitive(), geoIdJ.sensitive()}, sum});
      }
    }
    std::size_t nSourceLinks = sourceLinkContainer.size();

    Acts::ActsDynamicMatrix G =
        Acts::ActsDynamicMatrix::Zero(2 * nSourceLinks, 4);

    Eigen::VectorXd X(2 * nSourceLinks);
    double x0 = surfaceCenters.at(m_cfg.startSurfaceGeoId)(m_cfg.primaryIdx);
    for (std::size_t i = 0; i < nSourceLinks; i++) {
      const auto& ssl =
          sourceLinkContainer.at(i).template get<SimpleSourceLink>();
      const auto& geoId = ssl.geometryId();
      const auto* surface = m_cfg.surfaceMap.at(geoId);
      const Acts::Vector3& parameters =
          surface->localToGlobal(gctx, ssl.parametersLoc(), dir);
      X(2 * i) = parameters(m_cfg.longIdx);
      X(2 * i + 1) = parameters(m_cfg.shortIdx);

      G(2 * i, 0) = 1;
      G(2 * i, 2) = surface->center(gctx)(m_cfg.primaryIdx) - x0;
      G(2 * i + 1, 1) = 1;
      G(2 * i + 1, 3) = surface->center(gctx)(m_cfg.primaryIdx) - x0;
    }
    Eigen::LDLT<Acts::ActsDynamicMatrix> ldltD(constructCov(
        interSurfaceDistanceSums, sourceLinkContainer, dir, thetaMcsRms));
    Acts::SquareMatrix4 B = G.transpose() * ldltD.solve(G);
    Acts::Vector4 estimates = B.ldlt().solve(G.transpose() * ldltD.solve(X));

    pos(m_cfg.primaryIdx) = x0;
    pos(m_cfg.longIdx) = estimates(0);
    pos(m_cfg.shortIdx) = estimates(1);

    double tLong = estimates(2);
    double tShort = estimates(3);
    dir(m_cfg.primaryIdx) = 1.0;
    dir(m_cfg.longIdx) = tLong;
    dir(m_cfg.shortIdx) = tShort;
    dir.normalize();
    double denom = std::pow(1 + tLong * tLong + tShort * tShort, 1.5);

    Acts::ActsMatrix<eFreeSize, eEstimateSize> jacToGlob;
    jacToGlob.setZero();
    jacToGlob(eFreePos1, ePos0) = 1;
    jacToGlob(eFreePos2, ePos1) = 1;

    jacToGlob(eFreeDir0, eDir0) = -tLong / denom;
    jacToGlob(eFreeDir0, eDir1) = -tShort / denom;

    jacToGlob(eFreeDir1, eDir0) = (1 + tShort * tShort) / denom;
    jacToGlob(eFreeDir1, eDir1) = -tLong * tShort / denom;

    jacToGlob(eFreeDir2, eDir0) = -tLong * tShort / denom;
    jacToGlob(eFreeDir2, eDir1) = (1 + tLong * tLong) / denom;

    cov = jacToGlob * B.inverse() * jacToGlob.transpose();
    return (X - G * estimates).transpose() * ldltD.solve(X - G * estimates);
  }

  /// @brief get configuration instance
  const Config& config() const { return m_cfg; }

 private:
  /// @brief global measurement covariance construction
  ///
  /// @tparam container_t source link container type
  ///
  /// @param interSurfaceDistanceSums inter-surface distance dependent sums
  /// @param sourceLinksContainer source links container
  /// @param dir track direction estimate to be updated
  /// @param thetaMcsRms rms scattering angle in the layers
  ///
  /// @return global measurement covariance
  template <typename container_t>
  Acts::ActsDynamicMatrix constructCov(
      const std::unordered_map<std::pair<std::uint32_t, std::uint32_t>, double,
                               detail::PairHash>& interSurfaceDistanceSums,
      const container_t& sourceLinkContainer, const Acts::Vector3& dir,
      double thetaMcsRms) {
    std::size_t nSourceLinks = sourceLinkContainer.size();
    Acts::ActsDynamicMatrix D =
        Acts::ActsDynamicMatrix::Zero(2 * nSourceLinks, 2 * nSourceLinks);
    double mcsVarFactor = std::pow(dir(m_cfg.primaryIdx) * thetaMcsRms, 2);
    for (std::size_t i = 0; i < nSourceLinks; i++) {
      const auto& sslI =
          sourceLinkContainer.at(i).template get<SimpleSourceLink>();
      const auto& geoIdI = sslI.geometryId();

      double mcsVar =
          mcsVarFactor *
          interSurfaceDistanceSums.at({geoIdI.sensitive(), geoIdI.sensitive()});
      Acts::SquareMatrix2 varMat =
          sslI.covariance() + Acts::SquareMatrix2::Identity() * mcsVar;
      D.block(i * 2, i * 2, 2, 2) = varMat;
      for (std::size_t j = 0; j < i; j++) {
        const auto& sslJ =
            sourceLinkContainer.at(j).template get<SimpleSourceLink>();
        const auto& geoIdJ = sslJ.geometryId();
        double covPiPj =
            mcsVarFactor * interSurfaceDistanceSums.at(
                               {geoIdI.sensitive(), geoIdJ.sensitive()});

        Acts::SquareMatrix2 covMat = Acts::SquareMatrix2::Identity() * covPiPj;
        D.block(i * 2, j * 2, 2, 2) = covMat;
        D.block(j * 2, i * 2, 2, 2) = covMat;
      }
    }
    return std::move(D);
  }

  /// Configuration
  Config m_cfg;

  /// Measurement surfaces map
  std::vector<std::pair<Acts::GeometryIdentifier, const Acts::Surface*>>
      m_surfaces;
};
