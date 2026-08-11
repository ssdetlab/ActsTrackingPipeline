#include "TrackingPipeline/TrackFitting/KFTrackFittingAlgorithm.hpp"

#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <vector>

#include "TrackingPipeline/Utilities/NonOwningProjectionIterator.hpp"

KFTrackFittingAlgorithm::KFTrackFittingAlgorithm(const Config& config,
                                                 Acts::Logging::Level level)
    : IAlgorithm("TrackFittingAlgorithm", level), m_cfg(config) {
  m_inputTrackCandidates.initialize(m_cfg.inputTrackCandidates);
  m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);

  m_outputTrackContainer.initialize(m_cfg.outputTrackContainer);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ProcessCode KFTrackFittingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const auto& inputTrackCandidates = m_inputTrackCandidates(ctx);
  const auto& inputTrackParameters = m_inputTrackParameters(ctx);
  const auto& inputSourceLinks = m_inputSourceLinks(ctx);

  ACTS_DEBUG("Received " << inputTrackCandidates.size() << " track candidates");
  ACTS_DEBUG("Received " << inputTrackParameters.size() << " track parameters");
  ACTS_DEBUG("Received " << inputSourceLinks.size() << " source links");

  // Initialize KF options
  auto propOptions =
      KFFitterPropagatorOptions(ctx.geoContext, ctx.magFieldContext);
  propOptions.maxSteps = m_cfg.maxSteps;

  auto kfOptions = Acts::KalmanFitterOptions(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext, m_cfg.kfExtensions,
      propOptions);

  kfOptions.referenceSurface = m_cfg.referenceSurface;

  auto trackContainerBackend =
      std::make_shared<KFFitterTrackContainerBackend>();
  auto trackStateContainer = std::make_shared<KFFitterTrajectory>();
  Acts::TrackContainer trackContainer(trackContainerBackend,
                                      trackStateContainer);
  trackContainer.container().reserve(inputTrackCandidates.size());

  IndexTracks tracks;
  tracks.reserve(inputTrackCandidates.size());
  for (std::size_t idx = 0; idx < inputTrackCandidates.size(); idx++) {
    const auto& candidate = inputTrackCandidates.at(idx);
    const auto& sourceLinkIndices = candidate.sourceLinkIndices;
    if (sourceLinkIndices.empty()) {
      continue;
    }

    detail::NonOwningProjectionIterator begin(inputSourceLinks,
                                              sourceLinkIndices);
    detail::NonOwningProjectionIterator end = begin + sourceLinkIndices.size();
    const auto& startParameters =
        inputTrackParameters.at(candidate.originParametersIndex);

    auto res = m_cfg.fitter.fit(begin, end, startParameters, kfOptions,
                                trackContainer);
    if (res.ok()) {
      tracks.emplace_back(trackContainer.size() - 1,
                          candidate.originParametersIndex, candidate.trackId);
    }
  }
  tracks.shrink_to_fit();

  ACTS_DEBUG("Sending " << trackContainer.size() << " tracks");
  ACTS_DEBUG("Sending " << tracks.size() << " track index containers");

  m_outputTrackContainer(ctx, std::move(trackContainer));
  m_outputTracks(ctx, std::move(tracks));

  return ProcessCode::SUCCESS;
}
