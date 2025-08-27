#include "TrackingPipeline/TrackFinding/IdealSeedingAlgorithm.hpp"

#include "Acts/Utilities/Logger.hpp"

#include <cstddef>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

IdealSeedingAlgorithm::IdealSeedingAlgorithm(const Config& config,
                                             Acts::Logging::Level level)
    : IAlgorithm("IdealSeedingAlgorithm", level), m_cfg(config) {
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_inputSimClusters.initialize(m_cfg.inputSimClusters);
  m_outputSeeds.initialize(m_cfg.outputSeeds);
}

ProcessCode IdealSeedingAlgorithm::execute(const AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

  const std::vector<Acts::SourceLink>& inputSourceLinks =
      m_inputSourceLinks(ctx);
  const SimClusters& inputSimClusters = m_inputSimClusters(ctx);

  ACTS_DEBUG("Received " << inputSourceLinks.size() << " source links");

  if (inputSourceLinks.empty()) {
    ACTS_DEBUG("Input is empty. Skipping");
    m_outputSeeds(ctx, Seeds());
    return ProcessCode::SUCCESS;
  }

  using TrackId = std::tuple<int, int, int>;
  std::map<TrackId, std::vector<std::tuple<int, int, int>>> seedIdxs;

  Seeds outSeeds;
  for (std::size_t i = 0; i < inputSourceLinks.size(); i++) {
    const SimpleSourceLink& ssl =
        inputSourceLinks.at(i).get<SimpleSourceLink>();

    const SimCluster& cluster = inputSimClusters.at(ssl.index());
    for (std::size_t j = 0; j < cluster.truthHits.size(); j++) {
      const SimHit& hit = cluster.truthHits.at(j);
      seedIdxs[{hit.trackId, hit.parentTrackId, hit.runId}].push_back(
          {i, ssl.index(), j});
    }
  }
  ACTS_DEBUG("Found " << seedIdxs.size() << " seeds");

  for (const auto& [id, idxs] : seedIdxs) {
    if (idxs.size() < m_cfg.minSeedSize || idxs.size() > m_cfg.maxSeedSize) {
      continue;
    }

    std::set<Acts::GeometryIdentifier> seedGeoIds;
    std::vector<Acts::SourceLink> sourceLinks;
    sourceLinks.reserve(idxs.size());
    for (const auto& [slIdx, sslIdx, hitIdx] : idxs) {
      sourceLinks.push_back(inputSourceLinks.at(slIdx));
      seedGeoIds.insert(
          sourceLinks.back().get<SimpleSourceLink>().geometryId());
    }
    if (seedGeoIds.size() < m_cfg.minLayers ||
        seedGeoIds.size() > m_cfg.maxLayers) {
      continue;
    }

    outSeeds.emplace_back(std::move(sourceLinks),
                          inputSimClusters.at(std::get<1>(idxs.front()))
                              .truthHits.at(std::get<2>(idxs.front()))
                              .ipParameters,
                          static_cast<int>(outSeeds.size()));
  }
  ACTS_DEBUG("Sending " << outSeeds.size() << " seeds");
  m_outputSeeds(ctx, std::move(outSeeds));

  return ProcessCode::SUCCESS;
}
