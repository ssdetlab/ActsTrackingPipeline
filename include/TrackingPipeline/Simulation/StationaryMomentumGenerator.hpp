#pragma once

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"

/// @brief Stationary momentum generator
struct StationaryMomentumGenerator : public IMomentumGenerator {
  Acts::Vector3 momentum{0., 0., 0.};

  Acts::Vector3 genMomentum(RandomEngine& /*rng*/) const override {
    return momentum;
  }
};
