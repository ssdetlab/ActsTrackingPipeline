#pragma once

#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

/// @brief Uniform vertex generator
struct UniformVertexGenerator : public IVertexGenerator {
  Acts::Vector3 mins{0., 0., 0.};
  Acts::Vector3 maxs{0., 0., 0.};

  Acts::Vector3 genVertex(RandomEngine& rng) const override {
    std::uniform_real_distribution<double> uniform;
    Acts::Vector3 vertex{uniform(rng), uniform(rng), uniform(rng)};
    return mins + vertex.cwiseProduct(maxs - mins);
  }
};
