#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Interface for generating momentum vectors
struct IMomentumGenerator {
  virtual Acts::Vector3 genMomentum(RandomEngine& rng) const = 0;
};
