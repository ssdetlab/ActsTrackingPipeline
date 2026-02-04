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

#include <algorithm>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/TrackFitting/FittingServices.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmRegistry.hpp"

#include <toml.hpp>

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

double orthogonalLeastSquares(const std::vector<Acts::SourceLink>& sourceLinks,
                              Acts::Vector3& a, Acts::Vector3& b, int minGeoId,
                              int maxGeoId) {
  Acts::Vector3 meanVector(0, 0, 0);

  int nPoints = 0;
  for (std::size_t i = 0; i < sourceLinks.size(); i++) {
    const auto& sl = sourceLinks.at(i).get<SimpleSourceLink>();
    if (sl.geometryId().sensitive() < minGeoId ||
        sl.geometryId().sensitive() > maxGeoId) {
      continue;
    }
    nPoints++;
  }

  Eigen::MatrixXf points = Eigen::MatrixXf::Constant(nPoints, 3, 0);
  int k = 0;
  for (std::size_t i = 0; i < sourceLinks.size(); i++) {
    const auto& sl = sourceLinks.at(i).get<SimpleSourceLink>();
    if (sl.geometryId().sensitive() < minGeoId ||
        sl.geometryId().sensitive() > maxGeoId) {
      continue;
    }
    Acts::Vector3 parameters = sl.parametersGlob();
    points(k, 0) = parameters.x();
    points(k, 1) = parameters.y();
    points(k, 2) = parameters.z();

    meanVector += parameters;
    std::cout << "PARS " << parameters.transpose() << "\n";
    k++;
  }
  meanVector /= nPoints;
  a = meanVector;

  Eigen::MatrixXf centered = points.rowwise() - points.colwise().mean();
  Eigen::MatrixXf scatter = (centered.adjoint() * centered);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eig(scatter);
  Eigen::MatrixXf eigvecs = eig.eigenvectors();

  b[0] = eigvecs(0, 2);
  b[1] = eigvecs(1, 2);
  b[2] = eigvecs(2, 2);
  return eig.eigenvalues()(2);
}

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
    sourceLinkTrackContainer.push_back(candidate.sourceLinks);
    trackParametersContainer.push_back(candidate.ipParameters);
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

extern std::shared_ptr<AlignmentContext::AlignmentStore> g_alignmentStore;

namespace {

struct AlignmentAlgorithmRegistrar {
  AlignmentAlgorithmRegistrar() {
    using namespace TrackingPipeline;

    AlgorithmRegistry::instance().registerBuilder(
      "AlignmentAlgorithm",
      [](const toml::value& section,
         Acts::Logging::Level logLevel) -> AlgorithmPtr {

        auto& svc = FittingServices::instance();
        if (!svc.detector || !svc.baseKfOptions || !svc.referenceSurface || !svc.magneticField) {
          throw std::runtime_error(
              "AlignmentAlgorithm: FittingServices not initialized "
              "(detector/baseKfOptions/referenceSurface/magneticField missing)");
        }

        using Trajectory = FittingServices::Trajectory;
        Acts::KalmanFitterOptions<Trajectory> kfOptions = *svc.baseKfOptions;

        auto maxSteps = toml::find_or<unsigned int>(section, "maxSteps", 1000u);
        kfOptions.propagatorPlainOptions.maxSteps = maxSteps;
        kfOptions.referenceSurface = svc.referenceSurface.get();

        // Initial config
        AlignmentAlgorithm::Config cfg{
            /*inputTrackCandidates*/
            toml::find<std::string>(section, "inputTrackCandidates"),
            /*outputAlignmentParameters*/
            toml::find<std::string>(section, "outputAlignmentParameters"),
            /*align*/ nullptr,
            /*alignedTransformUpdater*/ ActsAlignment::AlignedTransformUpdater(),
            /*alignedDetElements*/ {},
            /*kfOptions*/ std::move(kfOptions),
            /*chi2ONdfCutOff*/ 1e-3,
            /*deltaChi2ONdfCutOff*/ {10, 1e-5},
            /*maxNumIterations*/ 200,
            /*alignmentMask*/
            ActsAlignment::AlignmentMask::Center1 |
            ActsAlignment::AlignmentMask::Center2 |
            ActsAlignment::AlignmentMask::Rotation2,
            /*alignmentMode*/ ActsAlignment::AlignmentMode::local};

        // Override chi2 / convergence from TOML if provided
        cfg.chi2ONdfCutOff =
            toml::find_or<double>(section, "chi2ONdfCutOff",
                                  cfg.chi2ONdfCutOff);
        std::size_t dChi2N =
            toml::find_or<std::size_t>(section, "deltaChi2ONdfCutOffN",
                                       cfg.deltaChi2ONdfCutOff.first);
        double dChi2Val =
            toml::find_or<double>(section, "deltaChi2ONdfCutOffValue",
                                  cfg.deltaChi2ONdfCutOff.second);
        cfg.deltaChi2ONdfCutOff = {dChi2N, dChi2Val};
        cfg.maxNumIterations =
            toml::find_or<std::size_t>(section, "maxNumIterations",
                                       cfg.maxNumIterations);

        // Optional mask/mode overrides
        const auto maskStr =
            toml::find_or<std::string>(section, "alignmentMask",
                                       "Center1Center2Rotation2");
        if (maskStr == "Center1Center2Rotation2") {
          cfg.alignmentMask = ActsAlignment::AlignmentMask::Center1 |
                              ActsAlignment::AlignmentMask::Center2 |
                              ActsAlignment::AlignmentMask::Rotation2;
        } else if (maskStr == "All") {
          cfg.alignmentMask = ActsAlignment::AlignmentMask::All;
        } else {
          throw std::runtime_error("Unknown alignmentMask '" + maskStr + "'");
        }

        const auto modeStr =
            toml::find_or<std::string>(section, "alignmentMode", "local");
        if (modeStr == "local") {
          cfg.alignmentMode = ActsAlignment::AlignmentMode::local;
        } else if (modeStr == "global") {
          cfg.alignmentMode = ActsAlignment::AlignmentMode::global;
        } else {
          throw std::runtime_error("Unknown alignmentMode '" + modeStr + "'");
        }

        // Alignment function
        auto detector = svc.detector;
        auto magneticField = svc.magneticField;
        cfg.align =
            AlignmentAlgorithm::makeAlignmentFunction(detector, magneticField);

        // all sensitive surfaces with id < 22 and != 10
        for (const auto& detElem : detector->detectorElements()) {
          auto* de = detElem.get();
          const auto& surf = de->surface();
          auto sid = surf.geometryId().sensitive();
          if (!sid) {
            continue;
          }
          if (sid < 22 && sid != 10) {
            cfg.alignedDetElements.push_back(de);
          }
        }

        // AlignedTransformUpdater: write transforms into the shared alignment store
        cfg.alignedTransformUpdater =
          [](Acts::DetectorElementBase* element,
             const Acts::GeometryContext&,
             const Acts::Transform3& trf) {
            if (!g_alignmentStore) {
              throw std::runtime_error(
                  "AlignmentAlgorithm: alignment store not set (g_alignmentStore is null)");
            }
            (*g_alignmentStore)[element->surface().geometryId()] = trf;
            return true;
          };

        return std::make_shared<AlignmentAlgorithm>(std::move(cfg), logLevel);
      });
  }
} _AlignmentAlgorithmRegistrar;

}  // namespace
