#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "TrackingPipeline/Utilities/ConcatVectorView.hpp"

namespace {

struct AlignmentFunctionImpl : public AlignmentAlgorithm::AlignmentFunction {
  Alignment align;

  explicit AlignmentFunctionImpl(Alignment&& a) : align(std::move(a)) {}

  AlignmentResult operator()(
      const AlignmentAlgorithm::SourceLinkContainer& sourceLinks,
      const AlignmentAlgorithm::TrackParametersContainer& initialParameters,
      const ActsAlignment::AlignmentOptions<KFFitterOptions>& options)
      override {
    return align.align(sourceLinks, initialParameters, options).value();
  }
};

}  // namespace

using namespace Acts::UnitLiterals;

std::shared_ptr<AlignmentAlgorithm::AlignmentFunction>
AlignmentAlgorithm::makeAlignmentFunction(
    const std::shared_ptr<const Acts::Experimental::Detector>& detector,
    const std::shared_ptr<const Acts::MagneticFieldProvider>& magneticField) {
  Stepper stepper(magneticField);
  Navigator::Config cfg;
  cfg.detector = detector.get();
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Navigator navigator(cfg, Acts::getDefaultLogger("AlignmentDetectorNavigator",
                                                  Acts::Logging::INFO));
  Propagator propagator(std::move(stepper), std::move(navigator));
  KFFitter trackFitter(
      std::move(propagator),
      Acts::getDefaultLogger("AlignmentKalmanFilter", Acts::Logging::INFO));
  Alignment alignment(std::move(trackFitter));

  // Build the alignment functions. owns the alignment object.
  return std::make_shared<AlignmentFunctionImpl>(std::move(alignment));
}

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

  m_outputAlignmentParameters.initialize(m_cfg.outputAlignmentParameters);
  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
}

ProcessCode AlignmentAlgorithm::execute(const AlgorithmContext& ctx) const {
  // NOTE: have to make a copy here to enable annealing
  // TODO: figure out a way to introduce annealing through contexts
  auto inputSourceLinks = m_inputSourceLinks(ctx);
  const auto& inputTrackCandidates = m_inputTrackCandidates(ctx);

  // NOTE: this copy is intentional -- stores updated track parameters
  auto inputTrackParameters = m_inputTrackParameters(ctx);

  if (inputTrackCandidates.empty()) {
    m_outputAlignmentParameters(ctx, {});
    m_outputTrackParameters(ctx, {});
    return ProcessCode::SUCCESS;
  }

  // Add constraints
  ConcatVectorView inputSourceLinksView(inputSourceLinks, m_cfg.constraints);

  std::size_t nCandidates = inputTrackCandidates.size();

  // Prepare the input track collection
  SourceLinkContainer sourceLinkContainer;
  sourceLinkContainer.reserve(nCandidates);

  std::vector<std::size_t> trackParametersIndicesContainer;
  trackParametersIndicesContainer.reserve(nCandidates);
  for (std::size_t i = 0; i < nCandidates; ++i) {
    // The list of hits and the initial start parameters
    const auto& candidate = inputTrackCandidates.at(i);

    sourceLinkContainer.emplace_back(inputSourceLinksView,
                                     candidate.sourceLinkIndices);
    trackParametersIndicesContainer.push_back(candidate.originParametersIndex);
  }
  TrackParametersContainer trackParametersContainer(
      inputTrackParameters, trackParametersIndicesContainer);

  // Set the alignment options
  ActsAlignment::AlignmentOptions<KFFitterOptions> alignOptions(
      m_cfg.kfOptions, m_cfg.alignmentTransformUpdater,
      m_cfg.alignmentParametersSolver, m_cfg.alignedDetElements,
      m_cfg.chi2ONdfCutOff, m_cfg.deltaChi2ONdfCutOff,
      m_cfg.maxAlignmentFitNumIt, m_cfg.alignmentMask);

  ACTS_DEBUG("Invoke track-based alignment with " << nCandidates
                                                  << " input tracks");

  AlignmentResult alignmentResult;
  for (std::size_t i = 0; i < m_cfg.nRefittingIt; i++) {
    if (m_cfg.annealingScheduler != nullptr) {
      double alpha = m_cfg.annealingScheduler->getAnnealingFactor(i);

      ACTS_INFO("Setting annealing factor "
                << alpha << " for alignment iteration " << i);

      // NOTE: workaround with the copy
      // TODO: find another way of applying annealing (contexts)
      for (auto& sl : inputSourceLinks) {
        auto& ssl = sl.get<SimpleSourceLink>();
        if (ssl.geometryId().sensitive() > 40) {
          ssl.setCovarianceAnnealingFactor(alpha);
        }
      }
    }

    ACTS_INFO("Staring alignment iteration " << i);
    auto result = (*m_cfg.alignmentFunction)(
        sourceLinkContainer, trackParametersContainer, alignOptions);
    ACTS_INFO("Alignment iteration " << i
                                     << ", deltaChi2 = " << result.deltaChi2);

    for (std::size_t j = 0; j < inputTrackCandidates.size(); j++) {
      if (m_cfg.annealingScheduler != nullptr) {
        // NOTE: workaround with the copy
        // TODO: find another way of applying annealing (contexts)
        for (auto& sl : inputSourceLinks) {
          auto& ssl = sl.get<SimpleSourceLink>();
          ssl.setCovarianceAnnealingFactor(1);
        }
      }

      const auto& candidate = inputTrackCandidates.at(j);
      auto& originParameters =
          inputTrackParameters.at(candidate.originParametersIndex);

      originParameters =
          m_cfg.trackParametersEstimator
              ->estimateParameters(ctx.geoContext, inputSourceLinks,
                                   candidate.sourceLinkIndices,
                                   originParameters.direction(),
                                   originParameters.position(ctx.geoContext))
              .value();
    }

    if (i == m_cfg.nRefittingIt - 1) {
      ACTS_INFO(
          "Alignment finished with chi2/ndf = " << result.averageChi2ONdf);
      alignmentResult = std::move(result);
    }
  }

  // Add alignment parameters to event store
  m_outputAlignmentParameters(ctx, std::move(alignmentResult));
  m_outputTrackParameters(ctx, std::move(inputTrackParameters));

  return ProcessCode::SUCCESS;
}
