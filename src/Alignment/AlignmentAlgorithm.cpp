#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"

#include "Acts/Detector/Detector.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/SourceLink.hpp>

#include <memory>
#include <utility>
#include <vector>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

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
      cfg, Acts::getDefaultLogger("DetectorNavigator", Acts::Logging::INFO));
  Propagator propagator(std::move(stepper), std::move(navigator));
  Fitter trackFitter(std::move(propagator));
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

  if (m_cfg.anchorSurface != nullptr) {
    SimpleSourceLink anchorSourceLink(
        Acts::Vector2::Zero(),
        m_cfg.anchorSurface->center(Acts::GeometryContext()),
        Acts::SquareMatrix2::Identity(), m_cfg.anchorSurface->geometryId(), -1,
        -1);
    m_anchorSourceLink = std::make_shared<Acts::SourceLink>(anchorSourceLink);
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
  for (std::size_t i = 0; i < numTracksUsed; ++i) {
    // The list of hits and the initial start parameters
    const auto& candidate = trackCandidates.at(i);
    sourceLinkTrackContainer.push_back(candidate.sourceLinks);
    if (m_cfg.anchorSurface != nullptr) {
      sourceLinkTrackContainer.back().push_back(*m_anchorSourceLink);
      trackParametersContainer.push_back(candidate.ipParameters);
    } else {
      trackParametersContainer.push_back(candidate.ipParameters);
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
