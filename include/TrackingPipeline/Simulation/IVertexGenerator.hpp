#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Interface for generating vertex positions
struct IVertexGenerator {
  virtual Acts::Vector3 genVertex(RandomEngine& rng) const = 0;
};
