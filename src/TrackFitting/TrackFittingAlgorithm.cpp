#include "TrackingPipeline/TrackFitting/TrackFittingAlgorithm.hpp"

ProcessCode TrackFittingAlgorithm::execute(const AlgorithmContext& ctx) const {
  // Get the input seeds
  // from the context
  auto inputCandidates = m_inputTrackCandidates(ctx);

  auto trackContainer = std::make_shared<TrackContainer>();
  auto trackStateContainer = std::make_shared<Trajectory>();
  Acts::TrackContainer tracks(trackContainer, trackStateContainer);

  for (const auto& candidate : inputCandidates) {
    auto start = candidate.ipParameters;
    auto sourceLinks = candidate.sourceLinks;

    auto res = m_cfg.fitter.fit(sourceLinks.begin(), sourceLinks.end(), start,
                                m_cfg.kfOptions, tracks);
  }
  m_outputTracks(ctx, std::move(tracks));

  return ProcessCode::SUCCESS;
}
