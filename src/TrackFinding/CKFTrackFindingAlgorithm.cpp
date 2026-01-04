#include "TrackingPipeline/TrackFinding/CKFTrackFindingAlgorithm.hpp"

#include <cstddef>
#include <memory>

// Note: Can send appropriate size candidates right away
// but will have to sacrifice the track view
ProcessCode CKFTrackFindingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Get the input seeds
  // from the context
  const auto& inputSeeds = m_inputSeeds(ctx);

  ACTS_DEBUG("Received " << inputSeeds.size() << " seeds");

  auto options = CKFOptions(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
      Acts::SourceLinkAccessorDelegate<SimpleSourceLinkAccessor::Iterator>{},
      m_cfg.extensions,
      Acts::PropagatorPlainOptions(ctx.geoContext, ctx.magFieldContext));
  options.propagatorPlainOptions.maxSteps = m_cfg.maxSteps;

  SimpleSourceLinkAccessor slAccessor;
  options.sourcelinkAccessor.template connect<&SimpleSourceLinkAccessor::range>(
      &slAccessor);

  auto container = std::make_shared<Acts::VectorTrackContainer>();
  auto trajectory = std::make_shared<Acts::VectorMultiTrajectory>();
  Acts::TrackContainer candidateContainer{container, trajectory};

  std::vector<int> trackIds;
  trackIds.reserve(inputSeeds.size());

  std::vector<Acts::CurvilinearTrackParameters> ipParametersGuesses;
  ipParametersGuesses.reserve(inputSeeds.size());

  Seeds trackCandidates;
  trackCandidates.reserve(inputSeeds.size());
  std::size_t idx = 0;
  for (const auto& seed : inputSeeds) {
    SimpleSourceLinkContainer ckfSourceLinks;
    for (auto& sl : seed.sourceLinks) {
      auto ssl = sl.get<SimpleSourceLink>();
      ckfSourceLinks.insert({ssl.geometryId(), std::move(ssl)});
    }
    slAccessor.container = &ckfSourceLinks;

    // Run the CKF for all initial track states
    const Acts::CurvilinearTrackParameters& ipParameters = seed.ipParameters;

    trackIds.push_back(seed.trackId);
    ipParametersGuesses.push_back(ipParameters);
    auto res = m_cfg.ckf.findTracks(ipParameters, options, candidateContainer);

    for (std::size_t tid = idx; tid < candidateContainer.size(); tid++) {
      const auto& track = candidateContainer.getTrack(tid);

      if (track.nOutliers() > 0) {
        continue;
      }
      std::vector<Acts::SourceLink> sourceLinks;
      for (const auto& trackState : track.trackStatesReversed()) {
        if (!trackState.hasUncalibratedSourceLink()) {
          continue;
        }
        sourceLinks.push_back(trackState.getUncalibratedSourceLink());
      }

      if (sourceLinks.size() < m_cfg.minCandidateSize ||
          sourceLinks.size() > m_cfg.maxCandidateSize) {
        continue;
      }

      trackCandidates.emplace_back(sourceLinks, ipParameters, tid);
    }
    idx = candidateContainer.size();
  }
  trackCandidates.shrink_to_fit();

  ACTS_DEBUG("Sending " << trackCandidates.size() << " track candidates");
  ACTS_DEBUG("Sending " << candidateContainer.size()
                        << " track candidates track views");
  m_outputTrackCandidates(ctx, std::move(trackCandidates));
  m_outputTrackView(ctx,
                    Tracks{std::move(candidateContainer), std::move(trackIds),
                           std::move(ipParametersGuesses)});

  return ProcessCode::SUCCESS;
}
