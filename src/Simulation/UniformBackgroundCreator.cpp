#include "TrackingPipeline/Simulation/UniformBackgroundCreator.hpp"

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

UniformBackgroundCreator::UniformBackgroundCreator(const Config& config)
    : m_cfg(config) {
  Acts::Vector2 stdDev = {m_cfg.resolution.first, m_cfg.resolution.second};
  m_cov = stdDev.cwiseProduct(stdDev).asDiagonal();
};

std::size_t UniformBackgroundCreator::gen(
    const AlgorithmContext& ctx, RandomEngine& rng, std::size_t id,
    std::vector<Acts::SourceLink>& sourceLinks, SimClusters& simClusters,
    std::vector<std::size_t>& sourceLinksIndices) const {
  int index = static_cast<int>(id);

  std::uniform_real_distribution uniform(0.0, 1.0);

  std::size_t currentSize = sourceLinksIndices.size();
  std::size_t resSize = m_cfg.surfaces.size() * m_cfg.nMeasurements;

  sourceLinks.reserve(currentSize + resSize);
  simClusters.reserve(currentSize + resSize);
  sourceLinksIndices.reserve(currentSize + resSize);

  for (const auto* surf : m_cfg.surfaces) {
    const auto& bounds = surf->bounds().values();

    for (std::size_t i = 0; i < m_cfg.nMeasurements; i++) {
      double x = bounds.at(0) + uniform(rng) * (bounds.at(2) - bounds.at(0));
      double y = bounds.at(1) + uniform(rng) * (bounds.at(3) - bounds.at(1));

      Acts::Vector2 hitLoc = Acts::Vector2(x, y);
      SimpleSourceLink ssl(
          hitLoc,
          surf->localToGlobal(ctx.geoContext, hitLoc, Acts::Vector3::UnitX()),
          m_cov, surf->geometryId(), ctx.eventNumber, sourceLinks.size());
      Acts::SourceLink sl(ssl);
      sourceLinksIndices.push_back(sourceLinks.size());
      sourceLinks.push_back(sl);

      SimCluster cluster{sl, {}, false};
      simClusters.push_back(cluster);
    }
  }
  return resSize;
};

const UniformBackgroundCreator::Config& UniformBackgroundCreator::config()
    const {
  return m_cfg;
}
