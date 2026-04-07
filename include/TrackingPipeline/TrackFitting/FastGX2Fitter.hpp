#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <cstddef>
#include <unordered_map>

/// @brief simple global chi2 fitter for intermediate estimations
///
/// The class performs a simple chi2 fit of a particle track
/// accounting for multiple scattering effects. The layers
/// of the tracking detector are assumed to be evenly spaced.
/// The fitter relies on the initial estimate of the track's direction.
class FastGX2Fitter {
 public:
  /// @brief config struct
  struct Config {
    /// Coordinate indices
    std::size_t primaryIdx;
    std::size_t longIdx;
    std::size_t shortIdx;
    /// Layer ids
    std::size_t firstLayerGeoId;
    std::size_t lastLayerGeoId;
    /// Surface map
    std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
        surfaceMap;
  };

  static const std::size_t covarianceDim = 6;
  using Covariance = Acts::ActsSquareMatrix<covarianceDim>;

  /// @brief constructor
  explicit FastGX2Fitter(const Config& cfg);

  /// @brief perform global chi2 fit
  ///
  /// @param gctx geometry context
  /// @param sourceLinks measurements participating in the fit
  /// @param pos track position estimate to be updated
  /// @param dir track direction estimate to be updated
  /// @param cov track covariance estimate to be updated
  ///
  /// @return chi2 of the track fit
  double gx2Fit(const Acts::GeometryContext& gctx,
                const std::vector<Acts::SourceLink>& sourceLinks,
                const std::vector<std::size_t>& sourceLinksIndices,
                double thetaMcsRms, Acts::Vector3& pos, Acts::Vector3& dir,
                Covariance& cov);

 private:
  /// @brief global covariance construction
  Acts::ActsDynamicMatrix constructCov(
      const std::vector<Acts::SourceLink>& sourceLinks,
      const std::vector<std::size_t>& sourceLinksIndices,
      const Acts::Vector3& dir, double thetaMcsRms);

  Config m_cfg;

  Acts::GeometryIdentifier m_firstLayerGeoId;
  Acts::GeometryIdentifier m_lastLayerGeoId;

  double m_primaryInterlayerDistance;
};
