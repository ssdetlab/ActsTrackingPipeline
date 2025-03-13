#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <random>
#include <stdexcept>
#include <vector>

#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

inline double powerLaw(double x) {
  const double scale = 5.31456e3;
  const double shift = -6.98928e1;
  const double pedestal = -4.99682e3;

  return (scale * std::pow(x + shift, -0.01) + pedestal) / 20011.3;
}

namespace E320Sim {

/// @brief Class that samples vertex from a ROOT histogram
struct E320PowerLawVertexGenerator : public IVertexGenerator {
  double yBoundLow;
  double yBoundHigh;

  double xBoundLow;
  double xBoundHigh;

  double yPower;
  double yShift;
  std::vector<double> zProbs;
  std::vector<double> zPositions;

  const double y0 = 90.3 - 29.94176 / 2;
  const double y1 = 331.1 + 29.94176 / 2;

  std::vector<double> is{y0, y1};
  std::vector<double> ws{powerLaw(y0), powerLaw(y1)};
    const double a = (powerLaw(y1) - powerLaw(y0)) / (y1 - y0);
    const double b = powerLaw(y0) - a * y0 + 1e-4;

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
