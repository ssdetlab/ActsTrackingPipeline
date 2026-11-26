#include "TrackingPipeline/TrackFitting/E320KFTrackFittingAlgorithm.hpp"

#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <optional>
#include <vector>

#include "TrackingPipeline/Utilities/NonOwningProjectionIterator.hpp"

E320::E320KFTrackFittingAlgorithm::E320KFTrackFittingAlgorithm(
    const Config& config, Acts::Logging::Level level)
    : IAlgorithm("TrackFittingAlgorithm", level), m_cfg(config) {
  m_inputTrackCandidates.initialize(m_cfg.inputTrackCandidates);
  m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);

  m_outputTrackContainer.initialize(m_cfg.outputTrackContainer);
  m_outputTracks.initialize(m_cfg.outputTracks);

  if (!m_cfg.mctxMagneticFields) {
    m_inputMagneticFieldParameters.initialize(
        m_cfg.inputMagneticFieldParameters);
  }
}

ProcessCode E320::E320KFTrackFittingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const auto& inputTrackCandidates = m_inputTrackCandidates(ctx);
  const auto& inputTrackParameters = m_inputTrackParameters(ctx);
  const auto& inputSourceLinks = m_inputSourceLinks(ctx);

  std::vector<std::shared_ptr<MagneticFieldStore>> dummyMagFieldContainer = {};
  const std::vector<std::shared_ptr<MagneticFieldStore>>&
      inputMagFieldParameters =
          m_cfg.mctxMagneticFields ? dummyMagFieldContainer
                                   : m_inputMagneticFieldParameters(ctx);

  ACTS_DEBUG("Received " << inputTrackCandidates.size() << " track candidates");
  ACTS_DEBUG("Received " << inputTrackParameters.size() << " track parameters");
  ACTS_DEBUG("Received " << inputSourceLinks.size() << " source links");

  // Initialize KF options
  auto propOptions =
      KFFitterPropagatorOptions(ctx.geoContext, ctx.magFieldContext);
  propOptions.maxSteps = m_cfg.maxSteps;
  propOptions.stepping.stepTolerance = m_cfg.stepTolerance;

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

  E320IndexTracks tracks;
  tracks.reserve(inputTrackCandidates.size());
  Acts::MagneticFieldContext candidateMctx{};
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

    if (!m_cfg.mctxMagneticFields) {
      candidateMctx =
          Acts::MagneticFieldContext{inputMagFieldParameters.at(idx)};
      kfOptions.magFieldContext =
          std::reference_wrapper<const Acts::MagneticFieldContext>(
              candidateMctx);
    }

    // Iterate the fit to linearize the propagation around the fitted q/p.
    // Only the q/p is updated, the position, direction and the prior
    // covariance are kept from the input track parameters
    Acts::CurvilinearTrackParameters iterParameters = startParameters;
    bool iterationsOk = true;
    for (std::size_t it = 1; it < m_cfg.nKFIterations; it++) {
      auto iterTrackContainerBackend =
          std::make_shared<KFFitterTrackContainerBackend>();
      auto iterTrackStateContainer = std::make_shared<KFFitterTrajectory>();
      Acts::TrackContainer iterTrackContainer(iterTrackContainerBackend,
                                              iterTrackStateContainer);

      auto iterRes = m_cfg.fitter.fit(begin, end, iterParameters, kfOptions,
                                      iterTrackContainer);
      if (!iterRes.ok()) {
        iterationsOk = false;
        break;
      }

      // Take the q/p of the smoothed state closest to the start of the fit
      std::optional<double> qOverP;
      for (const auto& state : iterRes.value().trackStatesReversed()) {
        if (state.hasSmoothed()) {
          qOverP = state.smoothed()[Acts::eBoundQOverP];
        }
      }
      if (!qOverP.has_value()) {
        iterationsOk = false;
        break;
      }
      ACTS_VERBOSE("Track " << candidate.trackId << " iteration " << it
                            << ": q/p " << iterParameters.qOverP() << " -> "
                            << qOverP.value());
      iterParameters.parameters()[Acts::eBoundQOverP] = qOverP.value();
    }
    if (!iterationsOk) {
      ACTS_DEBUG("Track " << candidate.trackId
                          << " failed in a KF iteration. Skip");
      continue;
    }

    auto res =
        m_cfg.fitter.fit(begin, end, iterParameters, kfOptions, trackContainer);
    if (res.ok()) {
      tracks.emplace_back(trackContainer.size() - 1,
                          candidate.originParametersIndex, idx,
                          candidate.trackId, candidate.htXCount);
    }
  }
  tracks.shrink_to_fit();

  ACTS_DEBUG("Sending " << trackContainer.size() << " tracks");
  ACTS_DEBUG("Sending " << tracks.size() << " track index containers");

  m_outputTrackContainer(ctx, std::move(trackContainer));
  m_outputTracks(ctx, std::move(tracks));

  return ProcessCode::SUCCESS;
}
