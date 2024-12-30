#pragma once

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

#include "Acts/Definitions/Algebra.hpp"

/// @brief Interface for generating momentum vectors
struct IMomentumGenerator {
    virtual Acts::Vector3 genMomentum(RandomEngine& rng) const = 0;
};
