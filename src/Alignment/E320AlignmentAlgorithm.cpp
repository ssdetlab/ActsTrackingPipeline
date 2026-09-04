#include "TrackingPipeline/Alignment/E320AlignmentAlgorithm.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

E320::E320AlignmentAlgorithm::E320AlignmentAlgorithm(const Config& cfg,
                                                     Acts::Logging::Level lvl)
    : IAlgorithm("E320::E320AlignmentAlgorithm", lvl), m_cfg(cfg) {
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

ProcessCode E320::E320AlignmentAlgorithm::execute(
    const AlgorithmContext& ctx) const {
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

  // Prepare indices containers
  std::vector<std::vector<std::size_t>> trackFitIndicesContainer;
  trackFitIndicesContainer.reserve(nCandidates);

  std::vector<std::vector<std::size_t>> initialTrackStateFitIndicesContainer;
  initialTrackStateFitIndicesContainer.reserve(nCandidates);

  std::vector<std::size_t> trackParametersIndicesContainer;
  trackParametersIndicesContainer.reserve(nCandidates);

  std::vector<std::size_t> magFieldIndicesContainer;
  magFieldIndicesContainer.reserve(nCandidates);
  for (std::size_t i = 0; i < nCandidates; ++i) {
    const auto& candidate = inputTrackCandidates.at(i);

    std::vector<std::size_t> trackFitIndices;
    trackFitIndices.reserve(candidate.sourceLinkIndices.size());

    std::vector<std::size_t> initialTrackStateFitIndices;
    initialTrackStateFitIndices.reserve(candidate.sourceLinkIndices.size());

    for (std::size_t idx : candidate.sourceLinkIndices) {
      const auto& geometryId =
          (inputSourceLinks.at(idx).type() == typeid(SimpleSourceLink))
              ? inputSourceLinks.at(idx).get<SimpleSourceLink>().geometryId()
              : inputSourceLinks.at(idx).get<ExtendedSourceLink>().geometryId();

      if (m_cfg.alignmentFitSurfaces.contains(geometryId)) {
        trackFitIndices.push_back(idx);
      }
      if (m_cfg.initialTrackStateFitSurfaces.contains(geometryId)) {
        initialTrackStateFitIndices.push_back(idx);
      }
    }

    trackFitIndicesContainer.push_back(std::move(trackFitIndices));
    initialTrackStateFitIndicesContainer.push_back(
        std::move(initialTrackStateFitIndices));
    trackParametersIndicesContainer.push_back(candidate.originParametersIndex);
    magFieldIndicesContainer.push_back(i);
  }

  // Initialize source link containers
  AlignmentFunction::SourceLinkContainer trackFitSourceLinkContainer;
  trackFitSourceLinkContainer.reserve(nCandidates);
  for (const auto& indices : trackFitIndicesContainer) {
    trackFitSourceLinkContainer.emplace_back(inputSourceLinks, indices);
  }

  AlignmentFunction::SourceLinkContainer
      initialTrackStateFitSourceLinkContainer;
  initialTrackStateFitSourceLinkContainer.reserve(nCandidates);
  for (const auto& indices : initialTrackStateFitIndicesContainer) {
    initialTrackStateFitSourceLinkContainer.emplace_back(inputSourceLinks,
                                                         indices);
  }

  AlignmentFunction::TrackParametersContainer trackParametersContainer(
      inputTrackParameters, trackParametersIndicesContainer);

  AlignmentFunction::MagneticFieldParametersContainer magFieldContainer(
      inputMagFieldParameters, magFieldIndicesContainer);

  // Run alignment
  ACTS_DEBUG("Invoke track-based alignment with " << nCandidates
                                                  << " input tracks");
  AlignmentFunction::AlignmentResult alignmentResult =
      (*m_cfg.alignmentFunction)(ctx.geoContext, ctx.magFieldContext,
                                 ctx.calibContext, trackFitSourceLinkContainer,
                                 initialTrackStateFitSourceLinkContainer,
                                 trackParametersContainer, magFieldContainer);
  ACTS_INFO(
      "Alignment finished with chi2/ndf = " << alignmentResult.averageChi2ONdf);

  // Collect updated track states
  auto outputTrackParameters = alignmentResult.updatedInitialTrackStates;

  // Add alignment parameters to event store
  m_outputAlignmentParameters(ctx, std::move(alignmentResult));
  m_outputTrackParameters(ctx, std::move(outputTrackParameters));

  return ProcessCode::SUCCESS;
}
