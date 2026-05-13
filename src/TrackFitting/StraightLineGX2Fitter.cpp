#include "TrackingPipeline/TrackFitting/StraightLineGX2Fitter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

#include <cstddef>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/detail/BinningValueUtils.hpp"

StraightLineGX2Fitter::StraightLineGX2Fitter(const Config& cfg) : m_cfg(cfg) {
  m_firstLayerGeoId.setSensitive(m_cfg.firstLayerGeoId);
  m_lastLayerGeoId.setSensitive(m_cfg.lastLayerGeoId);

  for (const auto& [geoId, surf] : m_cfg.surfaceMap) {
    m_surfaces.emplace_back(geoId, surf);
  }
  std::sort(m_surfaces.begin(), m_surfaces.end(),
            [](const auto& a, const auto& b) {
              return (a.first.sensitive() < b.first.sensitive());
            });
}

double StraightLineGX2Fitter::gx2Fit(
    const Acts::GeometryContext& gctx,
    const std::vector<Acts::SourceLink>& sourceLinks,
    const std::vector<std::size_t>& sourceLinksIndices, double thetaMcsRms,
    Acts::Vector3& pos, Acts::Vector3& dir, Acts::ActsSquareMatrix<6>& cov) {
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

  std::vector<std::size_t> filteredSourceLinksIndices;
  filteredSourceLinksIndices.reserve(m_cfg.surfaceMap.size());
  for (std::size_t i = 0; i < sourceLinksIndices.size(); i++) {
    std::size_t idx = sourceLinksIndices.at(i);
    const auto& ssl = sourceLinks.at(idx).get<SimpleSourceLink>();
    const auto& geoId = ssl.geometryId();
    if (geoId.sensitive() < m_firstLayerGeoId.sensitive() ||
        geoId.sensitive() > m_lastLayerGeoId.sensitive()) {
      continue;
    }
    filteredSourceLinksIndices.push_back(idx);
  }

  std::size_t nSourceLinks = filteredSourceLinksIndices.size();

  Acts::ActsDynamicMatrix G =
      Acts::ActsDynamicMatrix::Zero(2 * nSourceLinks, 4);

  Eigen::VectorXd X(2 * nSourceLinks);
  double x0 = surfaceCenters.at(m_firstLayerGeoId)(m_cfg.primaryIdx);
  for (std::size_t i = 0; i < nSourceLinks; i++) {
    std::size_t idx = filteredSourceLinksIndices.at(i);
    const auto& ssl = sourceLinks.at(idx).get<SimpleSourceLink>();
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
  Eigen::LDLT<Acts::ActsDynamicMatrix> ldltD(
      constructCov(interSurfaceDistanceSums, sourceLinks,
                   filteredSourceLinksIndices, dir, thetaMcsRms));
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

  Acts::ActsMatrix<6, 4> jacToGlob;
  jacToGlob.setZero();
  jacToGlob(1, 0) = 1;
  jacToGlob(2, 1) = 1;

  jacToGlob(3, 2) = -tLong / denom;
  jacToGlob(3, 3) = -tShort / denom;

  jacToGlob(4, 2) = (1 + tShort * tShort) / denom;
  jacToGlob(4, 3) = -tLong * tShort / denom;

  jacToGlob(5, 2) = -tLong * tShort / denom;
  jacToGlob(5, 3) = (1 + tLong * tLong) / denom;

  cov = jacToGlob * B.inverse() * jacToGlob.transpose();
  return (X - G * estimates).transpose() * ldltD.solve(X - G * estimates);
}

Acts::ActsDynamicMatrix StraightLineGX2Fitter::constructCov(
    const std::unordered_map<std::pair<std::uint32_t, std::uint32_t>, double,
                             detail::PairHash>& interSurfaceDistanceSums,
    const std::vector<Acts::SourceLink>& sourceLinks,
    const std::vector<std::size_t>& sourceLinksIndices,
    const Acts::Vector3& dir, double thetaMcsRms) {
  Acts::ActsDynamicMatrix D = Acts::ActsDynamicMatrix::Zero(
      2 * sourceLinksIndices.size(), 2 * sourceLinksIndices.size());
  double mcsVarFactor = std::pow(dir(m_cfg.primaryIdx) * thetaMcsRms, 2);
  for (std::size_t i = 0; i < sourceLinksIndices.size(); i++) {
    std::size_t idxI = sourceLinksIndices.at(i);
    const auto& sslI = sourceLinks.at(idxI).get<SimpleSourceLink>();
    const auto& geoIdI = sslI.geometryId();

    double mcsVar =
        mcsVarFactor *
        interSurfaceDistanceSums.at({geoIdI.sensitive(),
        geoIdI.sensitive()});
    Acts::SquareMatrix2 varMat =
        sslI.covariance() + Acts::SquareMatrix2::Identity() * mcsVar;
    D.block(i * 2, i * 2, 2, 2) = varMat;
    for (std::size_t j = 0; j < i; j++) {
      std::size_t idxJ = sourceLinksIndices.at(j);
      const auto& sslJ = sourceLinks.at(idxJ).get<SimpleSourceLink>();
      const auto& geoIdJ = sslJ.geometryId();
      double covPiPj =
          mcsVarFactor *
          interSurfaceDistanceSums.at({geoIdI.sensitive(),
          geoIdJ.sensitive()});

      Acts::SquareMatrix2 covMat = Acts::SquareMatrix2::Identity() * covPiPj;
      D.block(i * 2, j * 2, 2, 2) = covMat;
      D.block(j * 2, i * 2, 2, 2) = covMat;
    }
  }
  return std::move(D);
}
