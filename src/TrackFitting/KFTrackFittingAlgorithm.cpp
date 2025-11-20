#include "TrackingPipeline/TrackFitting/KFTrackFittingAlgorithm.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"

ProcessCode KFTrackFittingAlgorithm::execute(const AlgorithmContext& ctx) const {
  // Get the input seeds
  // from the context
  const auto& inputCandidates = m_inputTrackCandidates(ctx);

  ACTS_DEBUG("Received " << inputCandidates.size() << " track candidates");

  auto trackContainer = std::make_shared<TrackContainer>();
  auto trackStateContainer = std::make_shared<Trajectory>();
  Acts::TrackContainer tracks(trackContainer, trackStateContainer);

  std::vector<int> trackIds;
  trackIds.reserve(inputCandidates.size());

  std::vector<Acts::CurvilinearTrackParameters> ipParametersGuesses;
  ipParametersGuesses.reserve(inputCandidates.size());
  for (const auto& candidate : inputCandidates) {
    const auto& start = candidate.ipParameters;
    const auto& sourceLinks = candidate.sourceLinks;
    if (sourceLinks.empty()) {
      continue;
    }

    trackIds.push_back(candidate.trackId);
    ipParametersGuesses.push_back(candidate.ipParameters);

    auto res = m_cfg.fitter.fit(sourceLinks.begin(), sourceLinks.end(), start,
                                m_cfg.kfOptions, tracks);
  }
  trackIds.shrink_to_fit();
  ipParametersGuesses.shrink_to_fit();

  ACTS_DEBUG("Sending " << tracks.size() << " tracks");
  m_outputTracks(ctx, Tracks{tracks, trackIds, ipParametersGuesses});

  return ProcessCode::SUCCESS;
}
