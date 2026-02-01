#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Interface for generating momentum vectors
struct IMomentumGenerator {
  virtual Acts::Vector3 genMomentum(RandomEngine& rng) const = 0;

  virtual Acts::SquareMatrix4 getCovariance() const = 0;

  virtual Acts::Vector3 getMean() const = 0;
};
