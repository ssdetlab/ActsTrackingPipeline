#pragma once

#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/detail/NormalRandomVariable.hpp"

/// @brief Gaussian momentum generator
class GaussianVertexGenerator : public IVertexGenerator {
 public:
  struct Config {
    Acts::Vector3 mean;
    Acts::SquareMatrix3 cov;
  };

  GaussianVertexGenerator(const Config& cfg);

  Acts::Vector3 genVertex(RandomEngine& rng) const override;

  Acts::SquareMatrix3 getCovariance() const override;

  Acts::Vector3 getMean() const override;

 private:
  Config m_cfg;

  NormalRandomVariable m_normal;

  Acts::SquareMatrix3 m_cov;
};
