#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Interface for generating digitization smearing
struct IDigitizer {
  /// @brief digitize measurement
  ///
  /// @param rng random engine for sampling
  /// @param geoId geometry ID of the measurement surface
  /// @param pos 2D measurement
  ///
  /// @return digitized measurement
  virtual std::pair<Acts::Vector2, Acts::SquareMatrix2> digitize(
      RandomEngine& rng, const Acts::GeometryIdentifier& geoId,
      const Acts::Vector2& pos) const = 0;
};
