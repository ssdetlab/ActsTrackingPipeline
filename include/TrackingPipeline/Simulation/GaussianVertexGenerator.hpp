#pragma once

#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/detail/NormalRandomVariable.hpp"

/// @brief Gaussian momentum generator
struct GaussianVertexGenerator : public IVertexGenerator {
    NormalRandomVariable normal;

    GaussianVertexGenerator(Acts::Vector3 mean, Acts::SquareMatrix3 cov)
        : normal(mean, cov) {};

    Acts::Vector3 gen(RandomEngine& rng) const override {
        return normal.gen(rng);
    }
};
