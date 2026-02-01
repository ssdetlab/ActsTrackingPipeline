#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Interface for generating digitization smearing
struct IDigitizer {
  virtual Acts::Vector2 genCluster(RandomEngine& rng,
                                   const Acts::GeometryIdentifier& geoId,
                                   const Acts::Vector2& pos) const = 0;

  virtual Acts::SquareMatrix2 getCovariance(
      const Acts::GeometryIdentifier& geoId) const = 0;
};
