#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

/// @brief Uniform vertex generator
class UniformVertexGenerator : public IVertexGenerator {
 public:
  struct Config {
    Acts::Vector3 mins;
    Acts::Vector3 maxs;
  };

  UniformVertexGenerator(const Config& cfg);

  Acts::Vector3 genVertex(RandomEngine& rng) const override;

  Acts::SquareMatrix3 getCovariance() const override;

  Acts::Vector3 getMean() const override;

 private:
  Config m_cfg;

  Acts::SquareMatrix3 m_cov;

  Acts::Vector3 m_mean;
};
