#include "TrackingPipeline/TrackFinding/CKFTrackFindingAlgorithm.hpp"

ProcessCode CKFTrackFindingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Get the input seeds
  // from the context
  auto input = m_inputSeeds(ctx);

  auto options = CombinatorialKalmanFilterOptions(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
      Acts::SourceLinkAccessorDelegate<SimpleSourceLinkAccessor::Iterator>{},
      m_cfg.extensions,
      Acts::PropagatorPlainOptions(ctx.geoContext, ctx.magFieldContext));

  SimpleSourceLinkAccessor slAccessor;
  options.sourcelinkAccessor.template connect<&SimpleSourceLinkAccessor::range>(
      &slAccessor);

  TrackCandidates trackCandidates;
  for (const auto& seed : input) {
    SimpleSourceLinkContainer ckfSourceLinks;
    for (auto& sl : seed.sourceLinks) {
      auto ssl = sl.get<SimpleSourceLink>();
      ckfSourceLinks.insert({ssl.geometryId(), ssl});
    }

    slAccessor.container = &ckfSourceLinks;

    Acts::TrackContainer tc{Acts::VectorTrackContainer{},
                            Acts::VectorMultiTrajectory{}};

    // run the CKF for all initial track states
    Acts::CurvilinearTrackParameters ipParameters = seed.ipParameters;

    options.propagatorPlainOptions.maxSteps = 10000;

    auto res = m_cfg.ckf.findTracks(ipParameters, options, tc);
    if (!res.ok()) {
      continue;
    }

    for (std::size_t tid = 0u; tid < tc.size(); ++tid) {
      const auto track = tc.getTrack(tid);

      std::vector<Acts::SourceLink> sourceLinks;
      std::vector<double> predictedChi2;
      std::vector<double> filteredChi2;
      for (const auto trackState : track.trackStatesReversed()) {
        if (!trackState.hasUncalibratedSourceLink()) {
          continue;
        }
        sourceLinks.push_back(trackState.getUncalibratedSourceLink());
        predictedChi2.push_back(trackState.chi2());
        filteredChi2.push_back(trackState.chi2());
      }

      if (sourceLinks.size() < m_cfg.minCandidateSize ||
          sourceLinks.size() > m_cfg.maxCandidateSize) {
        continue;
      }

      trackCandidates.push_back(TrackCandidate{sourceLinks, ipParameters,
                                               seed.trackId, predictedChi2,
                                               filteredChi2});
    }
  }

  m_outputTrackCandidates(ctx, std::move(trackCandidates));

  return ProcessCode::SUCCESS;
}
