#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Interface for generating vertex positions
struct IVertexGenerator {
  virtual Acts::Vector3 genVertex(RandomEngine& rng) const = 0;

  virtual Acts::SquareMatrix3 getCovariance() const = 0;

  virtual Acts::Vector3 getMean() const = 0;
};
