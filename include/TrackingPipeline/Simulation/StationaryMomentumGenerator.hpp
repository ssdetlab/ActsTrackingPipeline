#pragma once

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"

/// @brief Stationary momentum generator
class StationaryMomentumGenerator : public IMomentumGenerator {
 public:
  struct Config {
    Acts::Vector3 momentum;
  };

  StationaryMomentumGenerator(const Config& config);

  Acts::Vector3 genMomentum(RandomEngine& rng) const override;

  Acts::SquareMatrix4 getCovariance() const override;

 private:
  Config m_cfg;

  Acts::SquareMatrix4 m_cov;
};
