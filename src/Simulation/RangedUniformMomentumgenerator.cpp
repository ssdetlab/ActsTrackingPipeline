#include "TrackingPipeline/Simulation/RangedUniformMomentumGenerator.hpp"

#include "Acts/Definitions/Algebra.hpp"

RangedUniformMomentumGenerator::RangedUniformMomentumGenerator(
    const Config& config)
    : m_cfg(config) {
  double meanP = 0;
  double mean2P = 0;
  for (const auto& [start, end] : m_cfg.pRanges) {
    meanP += (end + start) / 2.0;
    mean2P += (end * end + end * start + start * start) / 3.0;
  }
  meanP /= m_cfg.pRanges.size();
  mean2P /= m_cfg.pRanges.size();

  m_cov = Acts::SquareMatrix4::Zero();
  m_cov(3, 3) = (mean2P - meanP * meanP);
  m_mean = m_cfg.direction * meanP;

  m_rangeSelect =
      std::uniform_int_distribution<int>(0, m_cfg.pRanges.size() - 1);
}

Acts::Vector3 RangedUniformMomentumGenerator::genMomentum(
    RandomEngine& rng) const {
  int range = m_rangeSelect(rng);

  double pMin = m_cfg.pRanges.at(range).first;
  double pMax = m_cfg.pRanges.at(range).second;

  double p = pMin + m_uniform(rng) * (pMax - pMin);

  return p * m_cfg.direction;
}

Acts::SquareMatrix4 RangedUniformMomentumGenerator::getMomentumCovariance()
    const {
  return m_cov;
}

Acts::Vector3 RangedUniformMomentumGenerator::getMomentumMean() const {
  return m_mean;
}
