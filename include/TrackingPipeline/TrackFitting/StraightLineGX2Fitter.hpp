#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <cstddef>
#include <cstdint>
#include <unordered_map>

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
    /// Layer id range
    std::size_t firstLayerGeoId;
    std::size_t lastLayerGeoId;
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
  explicit StraightLineGX2Fitter(const Config& cfg);

  /// @brief perform global chi2 fit
  ///
  /// @param gctx geometry context
  /// @param sourceLinks source links container
  /// @param sourceLinks source links indices container
  /// @param thetaMcsRms rms scattering angle in the layers
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
  /// @brief global measurement covariance construction
  ///
  /// @param interSurfaceDistanceSums inter-surface distance dependent sums
  /// @param sourceLinks source links container
  /// @param sourceLinks source links indices container
  /// @param dir track direction estimate to be updated
  /// @param thetaMcsRms rms scattering angle in the layers
  ///
  /// @return global measurement covariance
  Acts::ActsDynamicMatrix constructCov(
      const std::unordered_map<std::pair<std::uint32_t, std::uint32_t>, double,
                               detail::PairHash>& interSurfaceDistanceSums,
      const std::vector<Acts::SourceLink>& sourceLinks,
      const std::vector<std::size_t>& sourceLinksIndices,
      const Acts::Vector3& dir, double thetaMcsRms);

  /// Configuration
  Config m_cfg;

  /// Measurement surfaces map
  std::vector<std::pair<Acts::GeometryIdentifier, const Acts::Surface*>>
      m_surfaces;

  /// Measurement surfaces geometry ID range
  Acts::GeometryIdentifier m_firstLayerGeoId;
  Acts::GeometryIdentifier m_lastLayerGeoId;
};
