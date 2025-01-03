#pragma once

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

#include "Acts/Definitions/Algebra.hpp"

/// @brief Interface for generating vertex positions
struct IVertexGenerator {
    virtual Acts::Vector3 genVertex(RandomEngine& rng) const = 0;
};
