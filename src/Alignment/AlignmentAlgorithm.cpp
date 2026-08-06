#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

AlignmentAlgorithm::AlignmentAlgorithm(const Config& cfg,
                                       Acts::Logging::Level lvl)
    : IAlgorithm("AlignmentAlgorithm", lvl), m_cfg(cfg) {
  if (m_cfg.inputTrackCandidates.empty()) {
    throw std::invalid_argument(
        "Missing input initial track parameters collection");
  }
  if (m_cfg.outputAlignmentParameters.empty()) {
    throw std::invalid_argument(
        "Missing output alignment parameters collection");
  }

  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_inputTrackCandidates.initialize(m_cfg.inputTrackCandidates);
  m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
  m_inputMagneticFieldParameters.initialize(m_cfg.inputMagneticFieldParameters);

  m_outputAlignmentParameters.initialize(m_cfg.outputAlignmentParameters);
  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
}

ProcessCode AlignmentAlgorithm::execute(const AlgorithmContext& ctx) const {
  const auto& inputSourceLinks = m_inputSourceLinks(ctx);
  const auto& inputTrackCandidates = m_inputTrackCandidates(ctx);
  const auto& inputTrackParameters = m_inputTrackParameters(ctx);
  const auto& inputMagFieldParameters = m_inputMagneticFieldParameters(ctx);

  if (inputTrackCandidates.empty()) {
    m_outputAlignmentParameters(ctx, {});
    m_outputTrackParameters(ctx, {});
    return ProcessCode::SUCCESS;
  }

  std::size_t nCandidates = inputTrackCandidates.size();

  // Prepare the input track collection
  SourceLinkContainer sourceLinkContainer;
  sourceLinkContainer.reserve(nCandidates);

  std::vector<std::size_t> trackParametersIndicesContainer;
  trackParametersIndicesContainer.reserve(nCandidates);
  for (std::size_t i = 0; i < nCandidates; ++i) {
    // The list of hits and the initial start parameters
    const auto& candidate = inputTrackCandidates.at(i);

    sourceLinkContainer.emplace_back(inputSourceLinks,
                                     candidate.sourceLinkIndices);
    trackParametersIndicesContainer.push_back(candidate.originParametersIndex);
  }
  TrackParametersContainer trackParametersContainer(
      inputTrackParameters, trackParametersIndicesContainer);

  ACTS_DEBUG("Invoke track-based alignment with " << nCandidates
                                                  << " input tracks");
  AlignmentResult alignmentResult = (*m_cfg.alignmentFunction)(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
      sourceLinkContainer, trackParametersContainer, inputMagFieldParameters);
  ACTS_INFO(
      "Alignment finished with chi2/ndf = " << alignmentResult.averageChi2ONdf);

  // TODO: remove
  auto outputTrackParameters = inputTrackParameters;

  // Add alignment parameters to event store
  m_outputAlignmentParameters(ctx, std::move(alignmentResult));
  m_outputTrackParameters(ctx, std::move(outputTrackParameters));

  return ProcessCode::SUCCESS;
}
