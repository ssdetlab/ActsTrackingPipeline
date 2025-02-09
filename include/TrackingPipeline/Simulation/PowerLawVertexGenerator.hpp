#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <random>
#include <stdexcept>
#include <vector>

#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

namespace E320Sim {

/// @brief Class that samples vertex from a ROOT histogram
struct E320PowerLawVertexGenerator : public IVertexGenerator {
  Acts::ActsScalar yBoundLow;
  Acts::ActsScalar yBoundHigh;

  Acts::ActsScalar xBoundLow;
  Acts::ActsScalar xBoundHigh;

  Acts::ActsScalar yPower;
  Acts::ActsScalar yShift;
  std::vector<Acts::ActsScalar> zProbs;
  std::vector<Acts::ActsScalar> zPositions;

  Acts::Vector3 genVertex(RandomEngine& rng) const override {
    if (zProbs.size() == 0) {
      throw std::runtime_error("Vertex generator not initialized");
    }

    std::discrete_distribution<> discrete(zProbs.begin(), zProbs.end());
    std::uniform_real_distribution<> uniform(0, 1);

    // Generate unfiorm x
    Acts::ActsScalar x = xBoundLow + (xBoundHigh - xBoundLow) * uniform(rng);

    // Generate power-law y
    Acts::ActsScalar gamma = uniform(rng);

    Acts::ActsScalar yTerm1 = gamma * std::pow(yBoundHigh - yShift, yPower + 1);
    Acts::ActsScalar yTerm2 =
        (1 - gamma) * std::pow(yBoundLow - yShift, yPower + 1);

    Acts::ActsScalar y = yShift + std::pow(yTerm1 + yTerm2, 1. / (yPower + 1));

    // Generate uniform z
    Acts::ActsScalar z = zPositions.at(discrete(rng));

    Acts::Vector3 glob(x, z, -y);

    return glob;
  }
};

}  // namespace E320Sim
