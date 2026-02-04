#include "TrackingPipeline/Simulation/GaussianVertexGenerator.hpp"

GaussianVertexGenerator::GaussianVertexGenerator(const Config& cfg)
    : m_cfg(cfg), m_normal(m_cfg.mean, m_cfg.cov) {};

Acts::Vector3 GaussianVertexGenerator::genVertex(RandomEngine& rng) const {
  return m_normal.gen(rng);
}

Acts::SquareMatrix3 GaussianVertexGenerator::getCovariance() const {
  return m_normal.getCovariance();
}

Acts::Vector3 GaussianVertexGenerator::getMean() const {
  return m_normal.getMean();
}
