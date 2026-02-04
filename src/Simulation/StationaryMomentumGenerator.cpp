#include "TrackingPipeline/Simulation/StationaryMomentumGenerator.hpp"

StationaryMomentumGenerator::StationaryMomentumGenerator(const Config& config)
    : m_cfg(config) {
  m_cov = Acts::SquareMatrix4::Zero();
}

Acts::Vector3 StationaryMomentumGenerator::genMomentum(
    RandomEngine& /*rng*/) const {
  return m_cfg.momentum;
}

Acts::SquareMatrix4 StationaryMomentumGenerator::getCovariance()
    const {
  return m_cov;
}
