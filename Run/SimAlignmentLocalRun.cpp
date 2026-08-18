#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Plugins/Json/JsonMaterialDecorator.hpp"
#include "Acts/Plugins/Json/MaterialMapJsonConverter.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <filesystem>
#include <iostream>
#include <memory>

#include <nlohmann/json.hpp>
#include <unistd.h>

#include "TrackingPipeline/Alignment/ActsAlignmentFunction.hpp"
#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"
#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Alignment/LocalAlignmentParametersSolverConstraints.hpp"
#include "TrackingPipeline/Alignment/LocalAlignmentParametersSolverSVD.hpp"
#include "TrackingPipeline/Alignment/LocalAlignmentTransformUpdater.hpp"
#include "TrackingPipeline/Alignment/detail/AlignmentStoreBuilders.hpp"
#include "TrackingPipeline/EventData/MixedSourceLinkCalibrator.hpp"
#include "TrackingPipeline/EventData/MixedSourceLinkSurfaceAccessor.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryOptions.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Infrastructure/TypeDefinitions.hpp"
#include "TrackingPipeline/Io/AlignmentParametersWriter.hpp"
#include "TrackingPipeline/Io/E320RootSimTrackReader.hpp"
#include "TrackingPipeline/Io/E320RootSimTrackWriter.hpp"
#include "TrackingPipeline/TrackFinding/E320TrackParametersEstimator.hpp"
#include "TrackingPipeline/TrackFitting/KFTrackFittingAlgorithm.hpp"
#include "toml++/toml.hpp"

using namespace Acts::UnitLiterals;

std::unique_ptr<const E320::GeometryOptions> E320::GeometryOptions::m_instance =
    nullptr;

int main() {
  // Geometry constraints instance
  const auto& goInst = *E320::GeometryOptions::instance();

  // Load configuration
  const std::string pathToCfg =
      "/home/romanurmanov/work/TrackingPipeline/ActsTrackingPipeline/conf/"
      "runs/"
      "SimAlignmentLocalRun.toml";
  auto runCfg = toml::parse_file(pathToCfg);

  auto getEntryDouble = [&runCfg](const std::string& section,
                                  const std::string& subsection) {
    return runCfg[section][subsection].value<double>().value();
  };
  auto getEntrySizeT = [&runCfg](const std::string& section,
                                 const std::string& subsection) {
    return runCfg[section][subsection].value<std::size_t>().value();
  };
  auto getEntryInt = [&runCfg](const std::string& section,
                               const std::string& subsection) {
    return runCfg[section][subsection].value<int>().value();
  };
  auto getEntryBool = [&runCfg](const std::string& section,
                                const std::string& subsection) {
    return runCfg[section][subsection].value<bool>().value();
  };
  auto getEntryStr = [&runCfg](const std::string& section,
                               const std::string& subsection) {
    return runCfg[section][subsection].value<std::string>().value();
  };

  // Set the log level
  Acts::Logging::Level logLevel =
      Acts::Logging::Level(getEntrySizeT("General", "logLevel"));

  // Initialize contexts
  Acts::GeometryContext gctx;

  // --------------------------------------------------------------
  // Detector setup

  // Material decorator
  Acts::MaterialMapJsonConverter::Config jsonMaterialConverterCfg;
  jsonMaterialConverterCfg.context = gctx;
  jsonMaterialConverterCfg.processSensitives =
      getEntryBool("MaterialMapJsonConverter", "processSensitives");
  jsonMaterialConverterCfg.processApproaches =
      getEntryBool("MaterialMapJsonConverter", "processApproaches");
  jsonMaterialConverterCfg.processRepresenting =
      getEntryBool("MaterialMapJsonConverter", "processRepresenting");
  jsonMaterialConverterCfg.processBoundaries =
      getEntryBool("MaterialMapJsonConverter", "processBoundaries");
  jsonMaterialConverterCfg.processVolumes =
      getEntryBool("MaterialMapJsonConverter", "processVolumes");
  jsonMaterialConverterCfg.processDenseVolumes =
      getEntryBool("MaterialMapJsonConverter", "processDenseVolumes");
  jsonMaterialConverterCfg.processNonMaterial =
      getEntryBool("MaterialMapJsonConverter", "processNonMaterial");

  auto materialDecorator = std::make_shared<Acts::JsonMaterialDecorator>(
      jsonMaterialConverterCfg,
      getEntryStr("JsonMaterialDecorator", "materialPath"), logLevel);

  // Construct detector
  auto detector = getEntryBool("Geometry", "materialDecorator")
                      ? E320::buildDetector(gctx, materialDecorator)
                      : E320::buildDetector(gctx, nullptr);

  // Set up surface maps for later use
  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
      gx2FitterSurfaceMap;
  std::unordered_set<Acts::GeometryIdentifier> alignmentFitSurfaces;
  std::unordered_set<Acts::GeometryIdentifier> initialTrackStateFitSurfaces;
  for (const auto& vol : detector->volumes()) {
    for (const auto& surf : vol->surfaces()) {
      const auto& geoId = surf->geometryId();
      if (geoId.sensitive() != 0u) {
        if (geoId.sensitive() >= goInst.tcParameters.front().geoId &&
            geoId.sensitive() <= goInst.tcParameters.back().geoId) {
          gx2FitterSurfaceMap[geoId] = surf;
          initialTrackStateFitSurfaces.insert(geoId);
          alignmentFitSurfaces.insert(geoId);
        }
      }
    }
  }

  // Initialize alignment store
  auto aStore = detail::makeAlignmentStore(gctx, detector.get());

  Acts::Vector3 globalShifts(0,
                             getEntryDouble("Geometry", "globalShiftY") * 1_mm,
                             getEntryDouble("Geometry", "globalShiftZ") * 1_mm);
  std::unordered_map<int, Acts::Vector3> localShifts{
      {10,
       Acts::Vector3(0_mm, getEntryDouble("Geometry", "localShiftY10") * 1_um,
                     getEntryDouble("Geometry", "localShiftZ10") * 1_um)},
      {12,
       Acts::Vector3(0_mm, getEntryDouble("Geometry", "localShiftY12") * 1_um,
                     getEntryDouble("Geometry", "localShiftZ12") * 1_um)},
      {14,
       Acts::Vector3(0_mm, getEntryDouble("Geometry", "localShiftY14") * 1_um,
                     getEntryDouble("Geometry", "localShiftZ14") * 1_um)},
      {16,
       Acts::Vector3(0_mm, getEntryDouble("Geometry", "localShiftY16") * 1_um,
                     getEntryDouble("Geometry", "localShiftZ16") * 1_um)},
      {18,
       Acts::Vector3(0_mm, getEntryDouble("Geometry", "localShiftY18") * 1_um,
                     getEntryDouble("Geometry", "localShiftZ18") * 1_um)}};
  Acts::Vector3 globalAngles(
      0_rad, 0_rad, getEntryDouble("Geometry", "globalAngleZ") * 1_mrad);
  std::unordered_map<int, Acts::Vector3> localAngles{
      {10, Acts::Vector3(0_rad, 0_rad,
                         getEntryDouble("Geometry", "localAngleZ10") * 1_mrad)},
      {12, Acts::Vector3(0_rad, 0_rad,
                         getEntryDouble("Geometry", "localAngleZ12") * 1_mrad)},
      {14, Acts::Vector3(0_rad, 0_rad,
                         getEntryDouble("Geometry", "localAngleZ14") * 1_mrad)},
      {16, Acts::Vector3(0_rad, 0_rad,
                         getEntryDouble("Geometry", "localAngleZ16") * 1_mrad)},
      {18,
       Acts::Vector3(0_rad, 0_rad,
                     getEntryDouble("Geometry", "localAngleZ18") * 1_mrad)}};
  aStore = detail::makeAlignmentStore(gctx, detector.get(), globalShifts,
                                      localShifts, globalAngles, localAngles);

  // Initialize alignment context
  AlignmentContext alignCtx(aStore);

  // Print alignment parameters
  Acts::GeometryContext defaultGctx;
  Acts::GeometryContext testCtx{alignCtx};
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive() != 0u) {
        std::cout << "-----------------------------------\n";
        std::cout << "SURFACE " << s->geometryId() << "\n";
        std::cout << "CENTER " << s->center(testCtx).transpose() << " -- "
                  << s->center(defaultGctx).transpose() << "\n";
        std::cout << "DELTA "
                  << (s->center(testCtx) - s->center(defaultGctx)).transpose() *
                         1e3
                  << "\n";
        std::cout << "NORMAL "
                  << s->normal(testCtx, s->center(testCtx),
                               Acts::Vector3::UnitY())
                         .transpose()
                  << " -- "
                  << s->normal(testCtx, s->center(defaultGctx),
                               Acts::Vector3::UnitY())
                         .transpose()
                  << "\n";
        std::cout << "ROTATION \n"
                  << s->transform(testCtx).rotation() << " -- \n"
                  << "\n"
                  << s->transform(defaultGctx).rotation() << "\n";

        std::cout << "EXTENT "
                  << s->polyhedronRepresentation(testCtx, 1000).extent()
                  << "\n -- \n"
                  << s->polyhedronRepresentation(defaultGctx, 1000).extent()
                  << "\n";
      }
    }
  }
  // Add alignment context to the geometry context
  gctx = Acts::GeometryContext{alignCtx};

  // --------------------------------------------------------------
  // The magnetic field setup

  auto field = E320::buildMagField(gctx);

  // --------------------------------------------------------------
  // Event reading
  SimpleSourceLink::SurfaceAccessor surfaceAccessor{detector.get()};

  // Setup the sequencer
  Sequencer::Config seqCfg;
  // seqCfg.events = 1e1;
  seqCfg.numThreads = 1;
  seqCfg.skip = 0;
  seqCfg.trackFpes = false;
  Sequencer sequencer(seqCfg);

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));

  // Add the sim data reader
  E320::E320RootSimTrackReader::Constraints readerConstraints{};
  readerConstraints.minSmoothedChi2 = 0;
  readerConstraints.maxSmoothedChi2 = 1e9;
  readerConstraints.minVertexEstLong = -std::numeric_limits<double>::max();
  readerConstraints.maxVertexEstLong = std::numeric_limits<double>::max();

  readerConstraints.minVertexEstShort = -std::numeric_limits<double>::max();
  readerConstraints.maxVertexEstShort = std::numeric_limits<double>::max();

  readerConstraints.minAbsMomentumEst = 0_GeV;
  readerConstraints.maxAbsMomentumEst = 10_GeV;

  E320::E320RootSimTrackReader::Config readerCfg;
  readerCfg.treeName = "fitted-tracks";
  readerCfg.outputSourceLinks = "SourceLinks";
  readerCfg.outputSimClusters = "SimClusters";
  readerCfg.outputSeedsGuess = "SeedsGuess";
  readerCfg.outputTrackParametersGuess = "TrackParametersGuess";
  readerCfg.outputSeedsEst = "SeedsEst";
  readerCfg.outputTrackParametersEst = "TrackParametersEst";
  readerCfg.outputMagneticFieldParameters = "MagFieldPars";
  readerCfg.constraints = readerConstraints;
  readerCfg.mergeIntoOneEvent = true;
  readerCfg.backwards = true;

  std::string pathToDir =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/sim/"
      "alignment/local/test";

  // Get the paths to the files in the directory
  for (const auto& entry : std::filesystem::directory_iterator(pathToDir)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".root") {
      continue;
    }
    std::string pathToFile = entry.path();
    readerCfg.filePaths.push_back(pathToFile);
  }

  // Add the reader to the sequencer
  sequencer.addReader(
      std::make_shared<E320::E320RootSimTrackReader>(readerCfg, logLevel));

  // --------------------------------------------------------------
  // Reference surface for sampling the track
  double halfX = std::numeric_limits<double>::max();
  double halfY = std::numeric_limits<double>::max();

  Acts::RotationMatrix3 refSurfToWorldRotationX =
      Acts::AngleAxis3(goInst.toWorldAngleX, Acts::Vector3::UnitX())
          .toRotationMatrix();
  Acts::RotationMatrix3 refSurfToWorldRotationY =
      Acts::AngleAxis3(goInst.toWorldAngleY, Acts::Vector3::UnitY())
          .toRotationMatrix();
  Acts::RotationMatrix3 refSurfToWorldRotationZ =
      Acts::AngleAxis3(goInst.toWorldAngleZ, Acts::Vector3::UnitZ())
          .toRotationMatrix();

  // Reestimation reference surface
  Acts::Transform3 reestimationRefSurfTransform = Acts::Transform3::Identity();
  reestimationRefSurfTransform.translation() = Acts::Vector3(
      goInst.ipTcDistance + 2 * goInst.tcHalfPrimary + 0.1_mm, 0, 0);
  reestimationRefSurfTransform.rotate(refSurfToWorldRotationX);
  reestimationRefSurfTransform.rotate(refSurfToWorldRotationY);
  reestimationRefSurfTransform.rotate(refSurfToWorldRotationZ);

  auto reestimationRefSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      reestimationRefSurfTransform,
      std::make_shared<Acts::RectangleBounds>(halfX, halfY));

  Acts::GeometryIdentifier reestimationRefSurfaceGeoId;
  reestimationRefSurfaceGeoId.setExtra(1);
  reestimationRefSurface->assignGeometryId(reestimationRefSurfaceGeoId);

  // Tracking reference surface
  Acts::Transform3 trackingRefSurfaceTransform = Acts::Transform3::Identity();
  trackingRefSurfaceTransform.translation() =
      Acts::Vector3(goInst.ipSurfaceCenterPrimary - 0.1_mm, 0, 0);
  trackingRefSurfaceTransform.rotate(refSurfToWorldRotationX);
  trackingRefSurfaceTransform.rotate(refSurfToWorldRotationY);
  trackingRefSurfaceTransform.rotate(refSurfToWorldRotationZ);

  auto trackingRefSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      trackingRefSurfaceTransform,
      std::make_shared<Acts::RectangleBounds>(halfX, halfY));

  Acts::GeometryIdentifier trackingRefSurfaceGeoId;
  trackingRefSurfaceGeoId.setExtra(2);
  trackingRefSurface->assignGeometryId(trackingRefSurfaceGeoId);

  // --------------------------------------------------------------
  // Alignment

  // Initialize track fitter extension
  KFFitterGainUpdater alignmentKFUpdater;
  KFFitterGainSmoother alignmentKFSmoother;

  Acts::KalmanFitterExtensions<KFFitterTrajectory> alignmentExtensions;
  // Add calibrator
  alignmentExtensions.calibrator
      .connect<&simpleSourceLinkCalibrator<KFFitterTrajectory>>();
  // Add the updater
  alignmentExtensions.updater
      .connect<&KFFitterGainUpdater::operator()<KFFitterTrajectory>>(
          &alignmentKFUpdater);
  // Add the smoother
  alignmentExtensions.smoother
      .connect<&KFFitterGainSmoother::operator()<KFFitterTrajectory>>(
          &alignmentKFSmoother);
  // Add the surface accessor
  alignmentExtensions.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  // Initial track state covariance matrix
  Acts::BoundVector trackOriginStdDevPrior;
  trackOriginStdDevPrior[Acts::eBoundLoc0] = 100_mm;
  trackOriginStdDevPrior[Acts::eBoundLoc1] = 100_mm;
  trackOriginStdDevPrior[Acts::eBoundPhi] = 10_degree;
  trackOriginStdDevPrior[Acts::eBoundTheta] = 10_degree;
  trackOriginStdDevPrior[Acts::eBoundQOverP] = 1 / 0.01_GeV;
  trackOriginStdDevPrior[Acts::eBoundTime] = 1_fs;
  Acts::BoundMatrix trackOriginCov =
      trackOriginStdDevPrior.cwiseProduct(trackOriginStdDevPrior).asDiagonal();

  // GX2 fitter setup
  Acts::GeometryIdentifier gx2StartSurfaceGeoId;
  gx2StartSurfaceGeoId.setSensitive(goInst.tcParameters.front().geoId);

  StraightLineGX2Fitter::Config gx2FitterCfg{};
  gx2FitterCfg.primaryIdx = goInst.primaryIdx;
  gx2FitterCfg.longIdx = goInst.longIdx;
  gx2FitterCfg.shortIdx = goInst.shortIdx;
  gx2FitterCfg.startSurfaceGeoId = gx2StartSurfaceGeoId;
  gx2FitterCfg.surfaceMap = gx2FitterSurfaceMap;

  auto gx2Fitter = std::make_shared<StraightLineGX2Fitter>(gx2FitterCfg);

  // Track parameters estimator
  E320::E320TrackParametersEstimator::Config trackParametersEstimatorCfg{};
  trackParametersEstimatorCfg.gx2Fitter = gx2Fitter;
  trackParametersEstimatorCfg.nIterations = 2;
  trackParametersEstimatorCfg.maxChi2 = std::numeric_limits<double>::max();
  trackParametersEstimatorCfg.referenceSurface = reestimationRefSurface.get();
  trackParametersEstimatorCfg.originCov =
      trackOriginStdDevPrior.cwiseProduct(trackOriginStdDevPrior).asDiagonal();
  trackParametersEstimatorCfg.propDirection =
      E320::E320TrackParametersEstimator::PropagationDirection::backward;

  auto trackParametersEstimator =
      std::make_shared<E320::E320TrackParametersEstimator>(
          trackParametersEstimatorCfg);

  // Alignment mask
  ActsAlignment::AlignmentMask alignmentMask =
      (ActsAlignment::AlignmentMask::Center1 |
       ActsAlignment::AlignmentMask::Center2 |
       ActsAlignment::AlignmentMask::Rotation2);

  // Alignment transform updater
  LocalAlignmentTransformUpdater::Config alignmentUpdaterCfg{};
  alignmentUpdaterCfg.alignmentStore = aStore;
  LocalAlignmentTransformUpdater alignmentUpdater(alignmentUpdaterCfg,
                                                  logLevel);

  // Alignment parameters solver

  LocalAlignmentParametersSolverConstraints::Config alignmentSolverCfg{};
  alignmentSolverCfg.alignmentMask = alignmentMask;
  LocalAlignmentParametersSolverConstraints alignmentSolver(alignmentSolverCfg,
                                                            logLevel);

  // LocalAlignmentParametersSolverSVD::Config alignmentSolverCfg{};
  // alignmentSolverCfg.alignmentMask = alignmentMask;
  // alignmentSolverCfg.maxSingularValueTol = 1e-5;
  // alignmentSolverCfg.singularValueGapTol = 9e-1;
  // LocalAlignmentParametersSolverSVD alignmentSolver(alignmentSolverCfg,
  //                                                   logLevel);

  // Number of refitting iterations
  std::size_t nRefittingIt = 2;

  // Alignment function
  ActsAlignmentFunction::Config alignmentFunctionCfg;
  alignmentFunctionCfg.detector = detector.get();
  alignmentFunctionCfg.magneticField = field;
  alignmentFunctionCfg.kfExtensions = alignmentExtensions;
  alignmentFunctionCfg.kfReferenceSurface = trackingRefSurface.get();
  alignmentFunctionCfg.maxKFSteps = static_cast<std::size_t>(1e5);
  alignmentFunctionCfg.chi2ONdfCutOff = 1e-16;
  alignmentFunctionCfg.deltaChi2ONdfCutOff = {10, 1e-5};
  alignmentFunctionCfg.maxAlignmentFitNumIt = 200;
  alignmentFunctionCfg.alignmentMask = alignmentMask;
  alignmentFunctionCfg.nRefittingIt = nRefittingIt;
  alignmentFunctionCfg.trackParametersEstimator = trackParametersEstimator;

  alignmentFunctionCfg.alignmentParametersSolver.connect<
      &LocalAlignmentParametersSolverConstraints::calculateAlignmentParameters>(
      &alignmentSolver);

  alignmentFunctionCfg.alignmentTransformUpdater
      .connect<&LocalAlignmentTransformUpdater::updateAlignmentParameters>(
          &alignmentUpdater);

  for (auto& det : detector->detectorElements()) {
    const auto& surface = det->surface();
    const auto& geoId = surface.geometryId().sensitive();
    if (geoId != 0u &&
        surface.geometryId().sensitive() > goInst.tcParameters.front().geoId &&
        surface.geometryId().sensitive() <= goInst.tcParameters.back().geoId) {
      alignmentFunctionCfg.alignedDetElements.push_back(det.get());
    }
  }

  auto alignmentFunction =
      std::make_shared<ActsAlignmentFunction>(alignmentFunctionCfg);

  // Alignment algorithm
  AlignmentAlgorithm::Config alignmentCfg;
  alignmentCfg.inputSourceLinks = "SourceLinks";
  alignmentCfg.inputTrackCandidates = "SeedsGuess";
  alignmentCfg.inputTrackParameters = "TrackParametersGuess";
  alignmentCfg.inputMagneticFieldParameters = "MagFieldPars";
  alignmentCfg.outputAlignmentParameters = "AlignmentParameters";
  alignmentCfg.outputTrackParameters = "UpdatedTrackParameters";
  alignmentCfg.alignmentFunction = alignmentFunction;
  alignmentCfg.alignmentFitSurfaces = alignmentFitSurfaces;
  alignmentCfg.initialTrackStateFitSurfaces = initialTrackStateFitSurfaces;

  auto alignmentAlgorithm =
      std::make_shared<AlignmentAlgorithm>(alignmentCfg, logLevel);
  sequencer.addAlgorithm(alignmentAlgorithm);

  // --------------------------------------------------------------
  // Track fitting
  KFFitterGainUpdater kfUpdater;
  KFFitterGainSmoother kfSmoother;

  // Initialize track fitter extensions
  Acts::KalmanFitterExtensions<KFFitterTrajectory> kfExtensions;
  // Add calibrator
  kfExtensions.calibrator
      .connect<&simpleSourceLinkCalibrator<KFFitterTrajectory>>();
  // Add the updater
  kfExtensions.updater
      .connect<&KFFitterGainUpdater::operator()<KFFitterTrajectory>>(
          &kfUpdater);
  // Add the smoother
  kfExtensions.smoother
      .connect<&KFFitterGainSmoother::operator()<KFFitterTrajectory>>(
          &kfSmoother);
  // Add the surface accessor
  kfExtensions.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  Navigator::Config cfg;
  cfg.detector = detector.get();
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Navigator kfNavigator(cfg,
                        Acts::getDefaultLogger("DetectorNavigator", logLevel));

  Acts::EigenStepper<> kfStepper(std::move(field));
  auto kfPropagator =
      Propagator(std::move(kfStepper), std::move(kfNavigator),
                 Acts::getDefaultLogger("Propagator", logLevel));

  const auto fitter = KFFitter(
      kfPropagator, Acts::getDefaultLogger("DetectorKalmanFilter", logLevel));

  // Add the track fitting algorithm to the sequencer
  KFTrackFittingAlgorithm::Config fitterCfg{
      .inputTrackCandidates = "SeedsGuess",
      .inputTrackParameters = "UpdatedTrackParameters",
      .inputSourceLinks = "SourceLinks",
      .outputTrackContainer = "TrackContainer",
      .outputTracks = "Tracks",
      .fitter = fitter,
      .maxSteps = static_cast<size_t>(1e5),
      .kfExtensions = kfExtensions,
      .referenceSurface = trackingRefSurface.get()};

  sequencer.addAlgorithm(
      std::make_shared<KFTrackFittingAlgorithm>(fitterCfg, logLevel));

  // --------------------------------------------------------------
  // Event write out

  // Fitted track writer
  E320::E320RootSimTrackWriter::Config trackWriterCfg;
  trackWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);
  trackWriterCfg.referenceSurface = trackingRefSurface.get();
  trackWriterCfg.inputTrackContainer = "TrackContainer";
  trackWriterCfg.inputTracks = "Tracks";
  trackWriterCfg.inputTrackParametersGuesses = "TrackParametersGuess";
  trackWriterCfg.inputSimClusters = "SimClusters";
  trackWriterCfg.treeName = "fitted-tracks";
  trackWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/sim/"
      "fitted-tracks.root";
  trackWriterCfg.applyPhiCorrection = true;

  sequencer.addWriter(
      std::make_shared<E320::E320RootSimTrackWriter>(trackWriterCfg, logLevel));

  // Alignment parameters writer
  AlignmentParametersWriter::Config alignmentWriterCfg;
  alignmentWriterCfg.treeName = "alignment-parameters";
  alignmentWriterCfg.inputAlignmentResults = "AlignmentParameters";
  alignmentWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/sim/"
      "alignment-parameters.root";

  sequencer.addWriter(std::make_shared<AlignmentParametersWriter>(
      alignmentWriterCfg, logLevel));

  sequencer.run();

  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive() != 0u) {
        std::cout << "-----------------------------------\n";
        std::cout << "SURFACE " << s->geometryId() << "\n";
        std::cout << "CENTER " << s->center(gctx).transpose() << " -- "
                  << s->center(Acts::GeometryContext()).transpose() << "\n";
        std::cout << "DELTA "
                  << (s->center(testCtx) - s->center(Acts::GeometryContext()))
                             .transpose() *
                         1e3
                  << "\n";

        std::cout << "NORMAL "
                  << s->normal(gctx, s->center(gctx), Acts::Vector3::UnitY())
                         .transpose()
                  << " -- "
                  << s->normal(Acts::GeometryContext(),
                               s->center(Acts::GeometryContext()),
                               Acts::Vector3::UnitY())
                         .transpose()
                  << "\n";
        std::cout << "ROTATION \n"
                  << s->transform(gctx).rotation() << " -- \n"
                  << "\n"
                  << s->transform(Acts::GeometryContext()).rotation() << "\n";
        std::cout << "EXTENT "
                  << s->polyhedronRepresentation(gctx, 1000).extent()
                  << "\n -- \n"
                  << s->polyhedronRepresentation(Acts::GeometryContext(), 1000)
                         .extent()
                  << "\n";
      }
    }
  }

  return 0;
}
