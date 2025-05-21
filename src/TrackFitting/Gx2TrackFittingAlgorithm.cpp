#include "Acts/EventData/TrackParameters.hpp"
#include <Acts/Utilities/Logger.hpp>

#include "TrackingPipeline/TrackFitting/Gx2TrackFittingAlgorithm.hpp"

ProcessCode Gx2TrackFittingAlgorithm::execute(const AlgorithmContext& ctx) const {
  // Get the input seeds
  // from the context
  auto inputCandidates = m_inputTrackCandidates(ctx);

  ACTS_DEBUG("Received " << inputCandidates.size() << " track candidates");

  auto trackContainer = std::make_shared<TrackContainer>();
  auto trackStateContainer = std::make_shared<Trajectory>();
  Acts::TrackContainer tracks(trackContainer, trackStateContainer);

  std::vector<int> trackIds;
  trackIds.reserve(inputCandidates.size());

  std::vector<Acts::CurvilinearTrackParameters> ipParametersGuesses;
  ipParametersGuesses.reserve(inputCandidates.size());
  for (const auto& candidate : inputCandidates) {
    auto start = candidate.ipParameters;
    auto sourceLinks = candidate.sourceLinks;

    trackIds.push_back(candidate.trackId);
    ipParametersGuesses.push_back(candidate.ipParameters);

    auto res = m_cfg.fitter.fit(sourceLinks.begin(), sourceLinks.end(), start,
                                m_cfg.options, tracks);
  }
  ACTS_DEBUG("Sending " << tracks.size() << " tracks");
  m_outputTracks(ctx, Tracks{tracks, trackIds, ipParametersGuesses});

  return ProcessCode::SUCCESS;
}
