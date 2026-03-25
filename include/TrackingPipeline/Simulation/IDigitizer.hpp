#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Interface for generating digitization smearing
struct IDigitizer {
  virtual std::pair<Acts::Vector2, Acts::SquareMatrix2> genCluster(
      RandomEngine& rng, const Acts::GeometryIdentifier& geoId,
      const Acts::Vector2& pos) const = 0;
};
