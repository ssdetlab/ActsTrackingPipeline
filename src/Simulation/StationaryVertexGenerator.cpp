#include "TrackingPipeline/Simulation/StationaryVertexGenerator.hpp"

StationaryVertexGenerator::StationaryVertexGenerator(const Config& cfg)
    : m_cfg(cfg) {
  m_cov = Acts::SquareMatrix3::Zero();
}

Acts::Vector3 StationaryVertexGenerator::genVertex(
    RandomEngine& /*rng*/) const {
  return m_cfg.vertex;
}

Acts::SquareMatrix3 StationaryVertexGenerator::getCovariance() const {
  return m_cov;
}

Acts::Vector3 StationaryVertexGenerator::getMean() const {
  return m_cfg.vertex;
}
