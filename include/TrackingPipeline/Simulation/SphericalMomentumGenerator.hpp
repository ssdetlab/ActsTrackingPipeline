#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"

/// @brief Spherical momentum generator
class SphericalMomentumGenerator : public IMomentumGenerator {
 public:
  struct Config {
    std::pair<double, double> pRange;
    std::pair<double, double> thetaRange;
    std::pair<double, double> phiRange;
  };

  SphericalMomentumGenerator(const Config& config);

  Acts::Vector3 genMomentum(RandomEngine& rng) const override;

  Acts::SquareMatrix4 getCovariance() const override;

  Acts::Vector3 getMean() const override;

 private:
  Config m_cfg;

  Acts::SquareMatrix4 m_cov;

  Acts::Vector3 m_mean;
};
