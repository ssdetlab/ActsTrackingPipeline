#include "TrackingPipeline/TrackFinding/TryAllSeedingAlgorithm.hpp"

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Utilities/IdxTree.hpp"

TryAllSeedingAlgorithm::TryAllSeedingAlgorithm(Config config,
                                               Acts::Logging::Level level)
    : IAlgorithm("TryAllTrackSeedingAlgorithm", level),
      m_cfg(std::move(config)) {
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_inputSourceLinkIndices.initialize(m_cfg.inputSourceLinkIndices);

  m_outputSeeds.initialize(m_cfg.outputSeeds);
  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
}

/// @brief The execute method
ProcessCode TryAllSeedingAlgorithm::execute(const AlgorithmContext& ctx) const {
  const auto& inputSourceLinks = m_inputSourceLinks(ctx);
  const auto& inputSourceLinkIndices = m_inputSourceLinkIndices(ctx);

  ACTS_DEBUG("Received " << inputSourceLinks.size() << " source links");
  ACTS_DEBUG("Received " << inputSourceLinkIndices.size()
                         << " source link indices");

  if (inputSourceLinkIndices.empty()) {
    ACTS_DEBUG("Input is empty. Skipping");
    m_outputSeeds(ctx, IndexSeeds());
    m_outputTrackParameters(ctx, {});
    return ProcessCode::SUCCESS;
  }

  IndexSeeds outSeeds;
  std::vector<Acts::CurvilinearTrackParameters> outTrackParameters;

  Acts::Vector3 point = Acts::Vector3::Zero();
  Acts::Vector3 dir = Acts::Vector3::UnitX();

  detail::IdxTree::IdxContainer idxContainer;
  idxContainer.reserve(inputSourceLinkIndices.size());
  std::size_t firstLayerId = std::numeric_limits<std::size_t>::max();
  for (std::size_t idx : inputSourceLinkIndices) {
    std::size_t geoId = inputSourceLinks.at(idx)
                            .get<SimpleSourceLink>()
                            .geometryId()
                            .sensitive();
    idxContainer.push_back({idx, geoId});
    if (firstLayerId > geoId) {
      firstLayerId = geoId;
    }
  }

  std::sort(idxContainer.begin(), idxContainer.end(),
            [](const auto& a, const auto& b) { return a.second < b.second; });
  auto rootEndIt = std::find_if(
      idxContainer.begin(), idxContainer.end(),
      [&firstLayerId](const auto& a) { return (a.second != firstLayerId); });

  for (auto it = idxContainer.begin(); it != rootEndIt; it++) {
    std::vector<std::size_t> trackContainer;
    trackContainer.reserve(m_cfg.minLayers);

    std::vector<std::vector<std::size_t>> splitSeedSlIdxs;
    detail::IdxTree idxTree(idxContainer, it, rootEndIt);
    idxTree.constructTracks(idxTree.getRootNode(), trackContainer,
                            splitSeedSlIdxs);
    for (const auto& seedIdxs : splitSeedSlIdxs) {
      if (seedIdxs.size() < m_cfg.minLayers ||
          seedIdxs.size() > m_cfg.maxLayers) {
        continue;
      }

      auto trackParametersResult =
          m_cfg.trackParametersEstimator->estimateParameters(
              ctx.geoContext, ctx.magFieldContext, inputSourceLinks, seedIdxs,
              dir, point);

      if (!trackParametersResult.ok()) {
        continue;
      }

      std::vector<std::size_t> sourceLinksIdxs = seedIdxs;
      sourceLinksIdxs.insert(sourceLinksIdxs.end(),
                             inputSourceLinkIndices.begin(),
                             inputSourceLinkIndices.end());

      outSeeds.emplace_back(std::move(sourceLinksIdxs),
                            outTrackParameters.size(), outSeeds.size());
      outTrackParameters.push_back(trackParametersResult.value());
    }
  }
  outSeeds.shrink_to_fit();
  outTrackParameters.shrink_to_fit();

  ACTS_DEBUG("Sending " << outSeeds.size() << " seeds");
  ACTS_DEBUG("Sending " << outTrackParameters.size() << " track parameteres");
  m_outputSeeds(ctx, std::move(outSeeds));
  m_outputTrackParameters(ctx, std::move(outTrackParameters));

  return ProcessCode::SUCCESS;
}

/// Get readonly access to the config parameters
const TryAllSeedingAlgorithm::Config& TryAllSeedingAlgorithm::config() const {
  return m_cfg;
}
