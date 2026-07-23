#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <unordered_map>
#include <utility>

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

/// @brief Surface-specific Gaussian-smearing digitizer
class SurfaceRangedDigitizer : public IDigitizer {
 public:
  /// Resolution shorthand
  using Resolution = std::pair<double, double>;

  /// @bried Nested config struct
  struct Config {
    /// Map of surface IDs to resolutions
    std::unordered_map<std::size_t, Resolution> resolutions;
  };

  /// @brief Constructor
  ///
  /// @param config configuration struct
  explicit SurfaceRangedDigitizer(const Config& config);

  /// @brief digitize measurement
  ///
  /// @param rng random engine for sampling
  /// @param geoId geometry ID of the measurement surface
  /// @param pos 2D measurement
  ///
  /// @return digitized measurement
  std::pair<Acts::Vector2, Acts::SquareMatrix2> digitize(
      RandomEngine& rng, const Acts::GeometryIdentifier& geoId,
      const Acts::Vector2& pos) const override;

 private:
  /// Configuration
  Config m_cfg;

  /// Measurement covariances
  std::unordered_map<std::size_t, Acts::SquareMatrix2> m_covs;

  /// Normal distribution instance
  mutable std::normal_distribution<double> m_normal{0., 1.};
};
