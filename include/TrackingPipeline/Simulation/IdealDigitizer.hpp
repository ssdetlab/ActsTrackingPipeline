#pragma once

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

#include "Acts/Definitions/Algebra.hpp"

struct IdealDigitizer : public IDigitizer {
    std::pair<Acts::SquareMatrix2, Acts::Vector2> gen(
        RandomEngine& rng,
        Acts::GeometryIdentifier /*geoId*/,
        Acts::Vector2 pos) const override {
            Acts::SquareMatrix2 cov = Acts::SquareMatrix2::Zero();
    
            return {cov, pos};
    }
};
