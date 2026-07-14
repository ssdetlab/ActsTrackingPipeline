#include "TrackingPipeline/Simulation/SurfaceRangedDigitizer.hpp"

#include <cstddef>

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
  std::size_t id = geoId.sensitive();

  Acts::Vector2 stdDev = {m_cfg.resolutions.at(id).first,
                          m_cfg.resolutions.at(id).second};
  Acts::Vector2 digLocal =
      pos + stdDev.cwiseProduct(Acts::Vector2(m_normal(rng), m_normal(rng)));
  return {digLocal, m_covs.at(id)};
}
