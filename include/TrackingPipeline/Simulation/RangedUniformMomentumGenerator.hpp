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

  explicit RangedUniformMomentumGenerator(const Config& config);

  Acts::Vector3 genMomentum(RandomEngine& rng) const override;

  Acts::SquareMatrix4 getMomentumCovariance() const override;

  Acts::Vector3 getMomentumMean() const override;

 private:
  Config m_cfg;

  Acts::SquareMatrix4 m_cov;

  Acts::Vector3 m_mean;

  mutable std::uniform_int_distribution<int> m_rangeSelect;
  mutable std::uniform_real_distribution<double> m_uniform{0, 1};
};
