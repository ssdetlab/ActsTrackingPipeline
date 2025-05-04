#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"

#include "Acts/EventData/TrackParameters.hpp"

AlignmentAlgorithm::AlignmentAlgorithm(Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("AlignmentAlgorithm", lvl), m_cfg(std::move(cfg)) {
  if (m_cfg.inputTrackCandidates.empty()) {
    throw std::invalid_argument(
        "Missing input initial track parameters collection");
  }
  if (m_cfg.outputAlignmentParameters.empty()) {
    throw std::invalid_argument(
        "Missing output alignment parameters collection");
  }

  m_inputTrackCandidates.initialize(m_cfg.inputTrackCandidates);
  m_outputAlignmentParameters.initialize(m_cfg.outputAlignmentParameters);
}

ProcessCode AlignmentAlgorithm::execute(const AlgorithmContext& ctx) const {
  // Read input data
  const auto& trackCandidates = m_inputTrackCandidates(ctx);

  std::size_t numTracksUsed = trackCandidates.size();

  // Prepare the input track collection
  std::vector<std::vector<Acts::SourceLink>> sourceLinkTrackContainer;
  sourceLinkTrackContainer.reserve(numTracksUsed);
  std::vector<Acts::CurvilinearTrackParameters> trackParametersContainer;
  trackParametersContainer.reserve(numTracksUsed);
  for (std::size_t itrack = 0; itrack < numTracksUsed; ++itrack) {
    // The list of hits and the initial start parameters
    const auto& candidate = trackCandidates.at(itrack);
    sourceLinkTrackContainer.push_back(candidate.sourceLinks);
    trackParametersContainer.push_back(candidate.ipParameters);
  }

  // Prepare the output for alignment parameters
  AlignmentParameters alignedParameters;

  // Set the alignment options
  ActsAlignment::AlignmentOptions<TrackFitterOptions> alignOptions(
      m_cfg.kfOptions, m_cfg.alignedTransformUpdater, m_cfg.alignedDetElements,
      m_cfg.chi2ONdfCutOff, m_cfg.deltaChi2ONdfCutOff, m_cfg.maxNumIterations);

  ACTS_DEBUG("Invoke track-based alignment with " << numTracksUsed
                                                  << " input tracks");
  auto result = (*m_cfg.align)(sourceLinkTrackContainer,
                               trackParametersContainer, alignOptions);
  if (result.ok()) {
    const auto& alignOutput = result.value();
    alignedParameters = alignOutput.alignedParameters;
    ACTS_VERBOSE(
        "Alignment finished with deltaChi2 = " << result.value().deltaChi2);
  } else {
    ACTS_WARNING("Alignment failed with " << result.error());
  }

  // add alignment parameters to event store
  m_outputAlignmentParameters(ctx, std::move(alignedParameters));
  return ProcessCode::SUCCESS;
}
