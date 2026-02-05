#include "TrackingPipeline/Simulation/UniformVertexGenerator.hpp"

#include "Acts/Definitions/Algebra.hpp"

UniformVertexGenerator::UniformVertexGenerator(const Config& cfg) : m_cfg(cfg) {
  m_cov = Acts::SquareMatrix3::Zero();
  m_cov(0, 0) = (m_cfg.maxs.x() - m_cfg.mins.x()) *
                (m_cfg.maxs.x() - m_cfg.mins.x()) / 12.0;
  m_cov(1, 1) = (m_cfg.maxs.y() - m_cfg.mins.y()) *
                (m_cfg.maxs.y() - m_cfg.mins.y()) / 12.0;
  m_cov(2, 2) = (m_cfg.maxs.z() - m_cfg.mins.z()) *
                (m_cfg.maxs.z() - m_cfg.mins.z()) / 12.0;

  m_mean(0) = (m_cfg.maxs.x() + m_cfg.mins.x()) / 2.0;
  m_mean(1) = (m_cfg.maxs.y() + m_cfg.mins.y()) / 2.0;
  m_mean(2) = (m_cfg.maxs.z() + m_cfg.mins.z()) / 2.0;
}

Acts::Vector3 UniformVertexGenerator::genVertex(RandomEngine& rng) const {
  std::uniform_real_distribution<double> uniform;
  Acts::Vector3 vertex{uniform(rng), uniform(rng), uniform(rng)};
  return m_cfg.mins + vertex.cwiseProduct(m_cfg.maxs - m_cfg.mins);
}

Acts::SquareMatrix3 UniformVertexGenerator::getCovariance() const {
  return m_cov;
}
Acts::Vector3 UniformVertexGenerator::getMean() const {
  return m_mean;
}
