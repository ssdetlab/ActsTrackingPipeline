#pragma once

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"

/// @brief Spherical momentum generator
struct SphericalMomentumGenerator : public IMomentumGenerator {
  std::pair<double, double> pRange;
  std::pair<double, double> thetaRange;
  std::pair<double, double> phiRange;

  Acts::Vector3 genMomentum(RandomEngine& rng) const override {
    std::uniform_real_distribution<double> uniform(0, 1);

    double pMag = pRange.first + (pRange.second - pRange.first) * uniform(rng);

    double phi =
        phiRange.first + (phiRange.second - phiRange.first) * uniform(rng);

    double theta =
        std::acos(std::cos(thetaRange.first) -
                  (std::cos(thetaRange.first) - std::cos(thetaRange.second)) *
                      uniform(rng));

    return pMag * Acts::Vector3(std::sin(theta) * std::cos(phi),
                                std::sin(theta) * std::sin(phi),
                                std::cos(theta));
  }
};
