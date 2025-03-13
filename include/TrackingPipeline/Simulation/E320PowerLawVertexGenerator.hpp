#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <random>
#include <stdexcept>
#include <vector>

#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

namespace E320Sim {

/// @brief Vertex sampler constructed to sample the inital
/// positions of the NCS-induced backgorund in the tracking
/// detector of the E320 experiment
struct E320PowerLawVertexGenerator : public IVertexGenerator {
  double yBoundLow;
  double yBoundHigh;

  double xBoundLow;
  double xBoundHigh;

  double yPower;
  double yShift;
  std::vector<double> zProbs;
  std::vector<double> zPositions;

  std::vector<double> is{yBoundLow, yBoundHigh};
  std::vector<double> ws{powerLaw(yBoundLow), powerLaw(yBoundHigh)};
  const double a =
      (powerLaw(yBoundHigh) - powerLaw(yBoundLow)) / (yBoundHigh - yBoundLow);
  const double b = powerLaw(yBoundLow) - a * yBoundLow + 1e-4;

  inline double powerLaw(double x) const {
    const double scale = 5.31456e3;
    const double shift = -6.98928e1;
    const double pedestal = -4.99682e3;

    return (scale * std::pow(x + shift, -0.01) + pedestal) / 20011.3;
  }
  Acts::Vector3 genVertex(RandomEngine& rng) const override {
    if (zProbs.empty() == 0) {
      throw std::runtime_error("Vertex generator not initialized");
    }

    std::discrete_distribution<> discrete(zProbs.begin(), zProbs.end());
    std::uniform_real_distribution<> uniform(0, 1);

    std::piecewise_linear_distribution<> linear(is.begin(), is.end(),
                                                ws.begin());

    // Generate unfiorm x
    double x = xBoundLow + (xBoundHigh - xBoundLow) * uniform(rng);

    // Generate power-law y
    bool sampled = false;
    double y;
    while (!sampled) {
      // Sample envelope
      double envY = linear(rng);

      // Sample ratio
      double ratio = uniform(rng);

      // Check ratio
      if (ratio < powerLaw(envY) / (a * envY + b)) {
        y = envY;
        sampled = true;
      }
    }

    // Generate uniform z
    double z = zPositions.at(discrete(rng));

    Acts::Vector3 glob(x, z, -y);

    return glob;
  }
};

}  // namespace E320Sim
