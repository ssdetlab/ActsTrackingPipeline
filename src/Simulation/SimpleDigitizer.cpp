#include "TrackingPipeline/Simulation/SimpleDigitizer.hpp"

SimpleDigitizer::SimpleDigitizer(const Config& config) : m_cfg(config) {
  Acts::Vector2 stdDev = {m_cfg.resolution.first, m_cfg.resolution.second};
  m_cov = stdDev.cwiseProduct(stdDev).asDiagonal();
}

Acts::Vector2 SimpleDigitizer::genCluster(
    RandomEngine& rng, const Acts::GeometryIdentifier& /*geoId*/,
    const Acts::Vector2& pos) const {
  std::normal_distribution<double> normal(0., 1.);

  Acts::Vector2 stdDev = {m_cfg.resolution.first, m_cfg.resolution.second};
  Acts::Vector2 digLocal =
      pos + stdDev.cwiseProduct(Acts::Vector2(normal(rng), normal(rng)));
  return digLocal;
}

Acts::SquareMatrix2 SimpleDigitizer::getCovariance(
    const Acts::GeometryIdentifier& /*geoId*/) const {
  return m_cov;
}
