#include "TrackingPipeline/Simulation/SurfaceRangedDigitizer.hpp"

SurfaceRangedDigitizer::SurfaceRangedDigitizer(const Config& config)
    : m_cfg(config) {
  for (const auto& [id, res] : m_cfg.resolutions) {
    Acts::Vector2 stdDev = {res.first, res.second};
    m_covs.insert({id, stdDev.cwiseProduct(stdDev).asDiagonal()});
  }
}

std::pair<Acts::Vector2, Acts::SquareMatrix2>
SurfaceRangedDigitizer::genCluster(RandomEngine& rng,
                                   const Acts::GeometryIdentifier& geoId,
                                   const Acts::Vector2& pos) const {
  std::normal_distribution<double> normal(0., 1.);

  Acts::Vector2 stdDev = {m_cfg.resolutions.at(geoId).first,
                          m_cfg.resolutions.at(geoId).second};
  Acts::Vector2 digLocal =
      pos + stdDev.cwiseProduct(Acts::Vector2(normal(rng), normal(rng)));
  return {digLocal, m_covs.at(geoId)};
}
