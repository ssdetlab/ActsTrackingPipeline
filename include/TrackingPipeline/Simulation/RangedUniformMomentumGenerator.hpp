#pragma once

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"

/// @brief Uniform momentum generator
struct RangedUniformMomentumGenerator : public IMomentumGenerator {
    std::vector<std::pair<Acts::ActsScalar, Acts::ActsScalar>> Pranges;

    Acts::Vector3 genMomentum(RandomEngine& rng) const override {
        std::uniform_int_distribution<int> range_select(0, Pranges.size() - 1);
        int range = range_select(rng);

        Acts::ActsScalar Pmin = Pranges.at(range).first;
        Acts::ActsScalar Pmax = Pranges.at(range).second;

        std::uniform_real_distribution<Acts::ActsScalar> uniform(Pmin, Pmax);
        Acts::ActsScalar p = uniform(rng);

        return p * Acts::Vector3(0, 1, 0);
    }
};
