#pragma once

#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

/// @brief Stationary vertex generator
struct StationaryVertexGenerator : public IVertexGenerator {
  Acts::Vector3 vertex{0., 0., 0.};

  Acts::Vector3 genVertex(RandomEngine& /*rng*/) const override {
    return vertex;
  }
};
