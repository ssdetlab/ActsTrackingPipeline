#include "TrackingPipeline/TrackFitting/TrackFittingAlgorithm.hpp"

ProcessCode TrackFittingAlgorithm::execute(const AlgorithmContext& ctx) const {
  // Get the input seeds
  // from the context
  auto input = m_inputSeeds(ctx);

  auto trackContainer = std::make_shared<TrackContainer>();
  auto trackStateContainer = std::make_shared<Trajectory>();
  Acts::TrackContainer tracks(trackContainer, trackStateContainer);

  std::vector<std::int32_t> trackIds;

  for (const auto& seed : input) {
    auto start = seed.ipParameters;
    auto sourceLinks = seed.sourceLinks;

    auto res = m_cfg.fitter.fit(sourceLinks.begin(), sourceLinks.end(), start,
                                m_cfg.kfOptions, tracks);

    if (!res.ok()) {
      ACTS_ERROR("Track fitting failed");
      continue;
    }

    trackIds.push_back(seed.trackId);
  }
  auto outTracks = Tracks<TrackContainer, Trajectory>{tracks, trackIds};

  m_outputTracks(ctx, std::move(outTracks));

  return ProcessCode::SUCCESS;
}
