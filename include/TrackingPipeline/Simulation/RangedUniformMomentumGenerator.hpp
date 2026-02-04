#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"

/// @brief Uniform momentum generator
class RangedUniformMomentumGenerator : public IMomentumGenerator {
 public:
  struct Config {
    std::vector<std::pair<double, double>> pRanges;
    Acts::Vector3 direction;
  };

  RangedUniformMomentumGenerator(const Config& config);

  Acts::Vector3 genMomentum(RandomEngine& rng) const override;

  Acts::SquareMatrix4 getCovariance() const override;

  Acts::Vector3 getMean() const override;

 private:
  Config m_cfg;

  Acts::SquareMatrix4 m_cov;

  Acts::Vector3 m_mean;
};
