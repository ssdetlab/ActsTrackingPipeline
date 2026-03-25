#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"

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

#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
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
    trackParametersContainer.push_back(candidate.ipParameters);
  }

  // Set the alignment options
  ActsAlignment::AlignmentOptions<TrackFitterOptions> alignOptions(
      m_cfg.kfOptions, m_cfg.alignedTransformUpdater, m_cfg.alignedDetElements,
      m_cfg.chi2ONdfCutOff, m_cfg.deltaChi2ONdfCutOff, m_cfg.maxNumIterations,
      m_cfg.alignmentMask, m_cfg.alignmentMode, m_cfg.maxSingularValueTol,
      m_cfg.singularValueGapTol, m_cfg.rigidAngleScale);

  ACTS_DEBUG("Invoke track-based alignment with " << numTracksUsed
                                                  << " input tracks");
  ActsAlignment::AlignmentResult alignmentResult;
  std::size_t nIt =
      (m_cfg.annealingScheduler == nullptr) ? 1 : m_cfg.nAnnealingIt;
  for (std::size_t i = 0; i < nIt; i++) {
    double alpha = (m_cfg.annealingScheduler == nullptr)
                       ? 1
                       : m_cfg.annealingScheduler->getAnnealingFactor(i);
    for (auto& sourceLinks : sourceLinkTrackContainer) {
      for (auto& sl : sourceLinks) {
        sl.get<SimpleSourceLink>().setCovarianceAnnealingFactor(alpha);
      }
    }

    ACTS_VERBOSE("Staring alignment iteration " << i << ": annealing factor "
                                                << alpha);
    auto result = (*m_cfg.align)(sourceLinkTrackContainer,
                                 trackParametersContainer, alignOptions);
    ACTS_VERBOSE("Alignment iteration "
                 << i << ": annealing factor " << alpha
                 << ", deltaChi2 = " << result.value().deltaChi2);
    if (i == nIt - 1) {
      if (result.ok()) {
        ACTS_VERBOSE(
            "Alignment finished with deltaChi2 = " << result.value().deltaChi2);
      } else {
        ACTS_WARNING("Alignment failed with " << result.error());
      }
      alignmentResult = std::move(result.value());
    }
  }

  m_outputAlignmentParameters(ctx, std::move(alignmentResult));

  // Add alignment parameters to event store
  return ProcessCode::SUCCESS;
}
