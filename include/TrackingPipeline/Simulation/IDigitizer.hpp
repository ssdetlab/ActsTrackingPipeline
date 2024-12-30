#pragma once

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

/// @brief Interface for generating digitization smearing
struct IDigitizer {
    virtual std::pair<Acts::SquareMatrix2, Acts::Vector2> genCluster(
        RandomEngine& rng,
        Acts::GeometryIdentifier geoId,
        Acts::Vector2 pos) const = 0;
};

