#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

/// @brief Stationary vertex generator
class StationaryVertexGenerator : public IVertexGenerator {
 public:
  struct Config {
    Acts::Vector3 vertex;
  };

  StationaryVertexGenerator(const Config& cfg);

  Acts::Vector3 genVertex(RandomEngine& rng) const override;

  Acts::SquareMatrix3 getCovariance() const override;

  Acts::Vector3 getMean() const override;

 private:
  Config m_cfg;

  Acts::SquareMatrix3 m_cov;
};
