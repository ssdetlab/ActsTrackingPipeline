#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

/// @brief Class digitizing measurements based
/// on the provided resolution assuming Gaussian smearing
class SimpleDigitizer : public IDigitizer {
 public:
  /// @bried Nested config struct
  struct Config {
    /// Digitization resolution
    std::pair<double, double> resolution;
  };

  /// @brief Constructor
  ///
  /// @param config configuration struct
  explicit SimpleDigitizer(const Config& config);

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

  /// Measurement covariance
  Acts::SquareMatrix2 m_cov;

  /// Normal distribution instance
  mutable std::normal_distribution<double> m_normal{0., 1.};
};
