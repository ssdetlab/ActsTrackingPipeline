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
}

Acts::Vector3 RangedUniformMomentumGenerator::genMomentum(
    RandomEngine& rng) const {
  std::uniform_int_distribution<int> range_select(0, m_cfg.pRanges.size() - 1);
  int range = range_select(rng);

  double Pmin = m_cfg.pRanges.at(range).first;
  double Pmax = m_cfg.pRanges.at(range).second;

  std::uniform_real_distribution<double> uniform(Pmin, Pmax);
  double p = uniform(rng);

  return p * m_cfg.direction;
}

Acts::SquareMatrix4 RangedUniformMomentumGenerator::getCovariance() const {
  return m_cov;
}

Acts::Vector3 RangedUniformMomentumGenerator::getMean() const {
  return m_mean;
}
