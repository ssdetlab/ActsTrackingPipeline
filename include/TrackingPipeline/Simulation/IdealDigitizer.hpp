#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

/// @brief Class digitizing measurements with infinite resolution
class IdealDigitizer : public IDigitizer {
 public:
  /// @brief Constructor
  IdealDigitizer() = default;

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
  /// Measurement covariance
  Acts::SquareMatrix2 m_cov;
};
