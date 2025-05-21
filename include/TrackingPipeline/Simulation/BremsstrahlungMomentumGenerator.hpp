#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cmath>
#include <random>

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/detail/NormalRandomVariable.hpp"

struct BremsstrahlungMomentumGenerator : public IMomentumGenerator {
  NormalRandomVariable normal;

  BremsstrahlungMomentumGenerator(Acts::Vector3 mean, Acts::SquareMatrix3 cov)
      : normal(mean, cov) {};

  Acts::Vector3 genMomentum(RandomEngine& rng) const override {
    Acts::Vector3 normalMom = normal.gen(rng);
    std::uniform_real_distribution<> uniform(0, 1);

    double Py = -std::log(1 - uniform(rng) * (1 - std::exp(-0.1 * 5))) / 0.1;
    return Acts::Vector3(normalMom.x(), Py, normalMom.z());
  }
};
