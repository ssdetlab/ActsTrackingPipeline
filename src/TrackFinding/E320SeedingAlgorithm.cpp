#include "TrackingPipeline/TrackFinding/E320SeedingAlgorithm.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <utility>
#include <vector>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"
#include "TrackingPipeline/Utilities/IdxTree.hpp"

namespace E320 {

E320SeedingAlgorithm::E320SeedingAlgorithm(const Config& config,
                                           Acts::Logging::Level level)
    : IAlgorithm("E320SeedingAlgorithm", level), m_cfg(config) {
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_inputDetSourceLinkIndices.initialize(m_cfg.inputDetSourceLinkIndices);
  m_inputBpmSourceLinkIndices.initialize(m_cfg.inputBpmSourceLinkIndices);

  m_outputSeeds.initialize(m_cfg.outputSeeds);
  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
}

ProcessCode E320SeedingAlgorithm::execute(const AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

  const auto& inputSourceLinks = m_inputSourceLinks(ctx);
  const auto& inputDetSourceLinkIndices = m_inputDetSourceLinkIndices(ctx);
  const auto& inputBpmSourceLinkIndices = m_inputBpmSourceLinkIndices(ctx);

  ACTS_DEBUG("Received " << inputSourceLinks.size() << " source links");
  ACTS_DEBUG("Received " << inputDetSourceLinkIndices.size()
                         << " detector source link indices");
  ACTS_DEBUG("Received " << inputBpmSourceLinkIndices.size()
                         << " bpm source link indices");

  if (inputDetSourceLinkIndices.empty()) {
    ACTS_DEBUG("Input is empty. Skipping");
    m_outputSeeds(ctx, IndexSeeds());
    m_outputTrackParameters(ctx, {});
    return ProcessCode::SUCCESS;
  }

  IndexSeeds outSeeds;
  std::vector<Acts::CurvilinearTrackParameters> outTrackParameters;
  std::vector<HoughTransformSeeder::HTSeed> htSeeds =
      m_cfg.htSeeder->findSeeds(ctx.geoContext, inputSourceLinks,
                                inputDetSourceLinkIndices, m_cfg.htOptions);

  ACTS_DEBUG("Found " << htSeeds.size() << " HT seeds");
  outSeeds.reserve(htSeeds.size() * m_cfg.maxLayers);
  outTrackParameters.reserve(htSeeds.size() * m_cfg.maxLayers);
  for (std::size_t i = 0; i < htSeeds.size(); i++) {
    const auto& [point, dir, slIdxs] = htSeeds.at(i);

    detail::IdxTree::IdxContainer idxContainer;
    idxContainer.reserve(slIdxs.size());
    std::size_t firstLayerId = std::numeric_limits<std::size_t>::max();
    for (std::size_t idx : slIdxs) {
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
                               inputBpmSourceLinkIndices.begin(),
                               inputBpmSourceLinkIndices.end());

        outSeeds.emplace_back(std::move(sourceLinksIdxs),
                              outTrackParameters.size(), outSeeds.size());
        outTrackParameters.push_back(trackParametersResult.value());
      }
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

}  // namespace E320
