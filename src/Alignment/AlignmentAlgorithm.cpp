#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

namespace {

using Updater = Acts::GainMatrixUpdater;
using Smoother = Acts::GainMatrixSmoother;
using Stepper = Acts::EigenStepper<>;
using Propagator =
    Acts::Propagator<Stepper, Acts::Experimental::DetectorNavigator>;
using Fitter = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;
using Alignment = ActsAlignment::Alignment<Fitter>;

struct AlignmentFunctionImpl : public AlignmentAlgorithm::AlignmentFunction {
  Alignment align;

  AlignmentFunctionImpl(Alignment&& a) : align(std::move(a)) {}

  AlignmentAlgorithm::AlignmentResult operator()(
      const std::vector<std::vector<Acts::SourceLink>>& sourceLinks,
      const std::vector<Acts::CurvilinearTrackParameters>& initialParameters,
      const ActsAlignment::AlignmentOptions<
          AlignmentAlgorithm::TrackFitterOptions>& options) override {
    return align.align(sourceLinks, initialParameters, options);
  };
};

}  // namespace

using namespace Acts::UnitLiterals;

std::shared_ptr<AlignmentAlgorithm::AlignmentFunction>
AlignmentAlgorithm::makeAlignmentFunction(
    const std::shared_ptr<const Acts::Experimental::Detector>& detector,
    const std::shared_ptr<const Acts::MagneticFieldProvider>& magneticField) {
  Stepper stepper(magneticField);
  Acts::Experimental::DetectorNavigator::Config cfg;
  cfg.detector = detector.get();
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Experimental::DetectorNavigator navigator(
      cfg, Acts::getDefaultLogger("AlignmentDetectorNavigator",
                                  Acts::Logging::INFO));
  Propagator propagator(std::move(stepper), std::move(navigator));
  Fitter trackFitter(
      std::move(propagator),
      Acts::getDefaultLogger("AlignmentKalmanFilter", Acts::Logging::INFO));
  Alignment alignment(std::move(trackFitter));

  // build the alignment functions. owns the alignment object.
  return std::make_shared<AlignmentFunctionImpl>(std::move(alignment));
}

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
  if (trackCandidates.empty()) {
    return ProcessCode::SUCCESS;
  }

  std::size_t numTracksUsed = trackCandidates.size();

  // Prepare the input track collection
  std::vector<std::vector<Acts::SourceLink>> sourceLinkTrackContainer;
  sourceLinkTrackContainer.reserve(numTracksUsed);
  std::vector<Acts::CurvilinearTrackParameters> trackParametersContainer;
  trackParametersContainer.reserve(numTracksUsed);
  for (std::size_t i = 0; i < numTracksUsed; ++i) {
    // The list of hits and the initial start parameters
    const auto& candidate = trackCandidates.at(i);

    std::vector<Acts::SourceLink> sourceLinks = candidate.sourceLinks;
    sourceLinks.insert(sourceLinks.end(), m_cfg.constraints.begin(),
                       m_cfg.constraints.end());
    sourceLinkTrackContainer.push_back(sourceLinks);

    if (m_cfg.propDirection == PropagationDirection::forward) {
      const auto& candidateIpPars = candidate.ipParameters;
      auto ipPars = Acts::CurvilinearTrackParameters(
          candidateIpPars.fourPosition(ctx.geoContext),
          candidateIpPars.direction(),
          candidateIpPars.charge() / candidateIpPars.absoluteMomentum(),
          m_cfg.originCov, candidateIpPars.particleHypothesis());
      trackParametersContainer.push_back(ipPars);
    } else {
      const auto& candidateIpPars = candidate.ipParameters;
      auto ipPars = Acts::CurvilinearTrackParameters(
          candidateIpPars.fourPosition(ctx.geoContext),
          -candidateIpPars.direction(),
          -candidateIpPars.charge() / candidateIpPars.absoluteMomentum(),
          m_cfg.originCov, candidateIpPars.particleHypothesis());
      trackParametersContainer.push_back(ipPars);
    }
  }

  // Prepare the output for alignment parameters
  AlignmentParameters alignedParameters;

  // Set the alignment options
  ActsAlignment::AlignmentOptions<TrackFitterOptions> alignOptions(
      m_cfg.kfOptions, m_cfg.alignedTransformUpdater, m_cfg.alignedDetElements,
      m_cfg.chi2ONdfCutOff, m_cfg.deltaChi2ONdfCutOff, m_cfg.maxNumIterations,
      m_cfg.alignmentMask, m_cfg.alignmentMode);

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

  // Add alignment parameters to event store
  m_outputAlignmentParameters(ctx, std::move(alignedParameters));
  return ProcessCode::SUCCESS;
}
