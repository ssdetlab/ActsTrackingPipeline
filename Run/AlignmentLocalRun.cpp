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

#include <Eigen/src/Core/ArithmeticSequence.h>
#include <nlohmann/json.hpp>
#include <unistd.h>

#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"
#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Alignment/LinearAnnealingScheduler.hpp"
#include "TrackingPipeline/Alignment/LocalAlignmentParametersSolverConstraints.hpp"
#include "TrackingPipeline/Alignment/LocalAlignmentParametersSolverSVD.hpp"
#include "TrackingPipeline/Alignment/LocalAlignmentTransformUpdater.hpp"
#include "TrackingPipeline/Alignment/detail/AlignmentStoreBuilders.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Infrastructure/TypeDefinitions.hpp"
#include "TrackingPipeline/Io/AlignmentParametersWriter.hpp"
#include "TrackingPipeline/Io/E320MagneticFieldParametersProvider.hpp"
#include "TrackingPipeline/Io/RootSeedWriter.hpp"
#include "TrackingPipeline/Io/RootTrackReader.hpp"
#include "TrackingPipeline/Io/RootTrackWriter.hpp"
#include "TrackingPipeline/MagneticField/MagneticFieldContextDecorator.hpp"
#include "TrackingPipeline/TrackFinding/E320TrackParametersEstimator.hpp"
#include "TrackingPipeline/TrackFitting/KFTrackFittingAlgorithm.hpp"

using namespace Acts::UnitLiterals;

std::unique_ptr<const E320::GeometryOptions> E320::GeometryOptions::m_instance =
    nullptr;

int main() {
  const auto& goInst = *E320::GeometryOptions::instance();

  // Set the log level
  Acts::Logging::Level logLevel = Acts::Logging::DEBUG;

  // Initialize contexts
  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext cctx;

  // --------------------------------------------------------------
  // Detector setup

  // Material decorator
  Acts::MaterialMapJsonConverter::Config jsonMaterialConverterCfg;
  jsonMaterialConverterCfg.context = gctx;
  jsonMaterialConverterCfg.processSensitives = true;
  jsonMaterialConverterCfg.processApproaches = true;
  jsonMaterialConverterCfg.processRepresenting = true;
  jsonMaterialConverterCfg.processBoundaries = true;
  jsonMaterialConverterCfg.processVolumes = true;
  jsonMaterialConverterCfg.processDenseVolumes = false;
  jsonMaterialConverterCfg.processNonMaterial = false;

  auto materialDecorator = std::make_shared<Acts::JsonMaterialDecorator>(
      jsonMaterialConverterCfg,
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_material/"
      "Uniform_DirectZ_Tracker_PDCWindow_256x128_1e6/material.json",
      logLevel);

  // Construct detector
  auto detector = E320::buildDetector(gctx, materialDecorator);

  // Set up surface maps for later use
  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
      gx2FitterSurfaceMap;
  for (const auto& vol : detector->volumes()) {
    std::cout << "------------------------------------------\n";
    std::cout << vol->name() << "\n";
    std::cout << vol->extent(gctx);
    std::cout << "Surfaces:\n";
    for (const auto& surf : vol->surfaces()) {
      std::cout << surf->geometryId() << "\n";
      std::cout << surf->center(gctx) << "\n";
      std::cout << surf->polyhedronRepresentation(gctx, 1000).extent() << "\n";
      if (surf->geometryId().sensitive() >= goInst.tcParameters.front().geoId &&
          surf->geometryId().sensitive() <= goInst.tcParameters.back().geoId) {
        gx2FitterSurfaceMap[surf->geometryId()] = surf;
      }
    }
  }

  // Initialize alignment store
  auto aStore = detail::makeAlignmentStore(gctx, detector.get());

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

  // E320::E320MagneticFieldParametersProvider::Config magFieldProviderCfg;
  // magFieldProviderCfg.treeName = "magnets";
  // std::vector<std::string> inDirs = {
  //     "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/data/"
  //     "alignment/local_march_2026_data/sig/magnets"};

  // // Get the paths to the files in the directory
  // for (const auto& dir : inDirs) {
  //   for (const auto& entry : std::filesystem::directory_iterator(dir)) {
  //     if (!entry.is_regular_file() || entry.path().extension() != ".root") {
  //       continue;
  //     }
  //     std::string pathToFile = entry.path();
  //     magFieldProviderCfg.filePaths.push_back(pathToFile);
  //   }
  // }
  // E320::E320MagneticFieldParametersProvider magFieldProvider(
  //     magFieldProviderCfg);

  // auto magFieldCollection =
  // magFieldProvider.getMagneticFieldStoreCollection();

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
  // sequencer.addContextDecorator(
  //     std::make_shared<MagneticFieldContextDecorator>(magFieldCollection));

  // Add the sim data reader
  RootTrackReader::Constraints readerConstraints{};
  readerConstraints.minChi2 = 0;
  readerConstraints.maxChi2 = 40;
  readerConstraints.minVertexEstLong = -2000;
  readerConstraints.maxVertexEstLong = 2000;

  readerConstraints.minVertexEstShort = -2000;
  readerConstraints.maxVertexEstShort = 2000;

  readerConstraints.minAbsMomentumEst = 0_GeV;
  readerConstraints.maxAbsMomentumEst = 10_GeV;

  RootTrackReader::Config readerCfg;
  readerCfg.treeName = "fitted-tracks";
  readerCfg.outputSourceLinks = "SourceLinks";
  readerCfg.outputSeedsGuess = "SeedsGuess";
  readerCfg.outputTrackParametersGuess = "TrackParametersGuess";
  readerCfg.outputSeedsEst = "SeedsEst";
  readerCfg.outputTrackParametersEst = "TrackParametersEst";
  readerCfg.constraints = readerConstraints;
  readerCfg.mergeIntoOneEvent = true;

  std::string pathToDir =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/data/"
      "alignment/local_march_2026_data/sig/5um_errors/"
      "it2/filtered";

  // Get the paths to the files in the directory
  for (const auto& entry : std::filesystem::directory_iterator(pathToDir)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".root") {
      continue;
    }
    std::string pathToFile = entry.path();
    readerCfg.filePaths.push_back(pathToFile);
  }

  // Add the reader to the sequencer
  sequencer.addReader(std::make_shared<RootTrackReader>(readerCfg, logLevel));

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
  reestimationRefSurfTransform.translation() =
      Acts::Vector3(goInst.ipTcDistance - 0.3_mm, 0, 0);
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
  trackingRefSurfaceTransform.translation() = Acts::Vector3(
      goInst.dipoleCenterPrimary + goInst.dipoleHalfPrimary + 0.01_mm, 0, 0);
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

  // Initialize track fitter extensions
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

  auto alignmentPropOptions = KFFitterPropagatorOptions(gctx, mctx);

  alignmentPropOptions.maxSteps = 1000;

  auto alignmentKFOptions = KFFitterOptions(
      gctx, mctx, cctx, alignmentExtensions, alignmentPropOptions);

  alignmentKFOptions.referenceSurface = trackingRefSurface.get();

  // Initial track state covariance matrix
  Acts::BoundVector trackOriginStdDevPrior;
  trackOriginStdDevPrior[Acts::eBoundLoc0] = 100_mm;
  trackOriginStdDevPrior[Acts::eBoundLoc1] = 100_mm;
  trackOriginStdDevPrior[Acts::eBoundTime] = 25_ns;
  trackOriginStdDevPrior[Acts::eBoundPhi] = 10_rad;
  trackOriginStdDevPrior[Acts::eBoundTheta] = 10_rad;
  trackOriginStdDevPrior[Acts::eBoundQOverP] = 1 / 0.01_GeV;
  Acts::BoundMatrix trackOriginCov =
      trackOriginStdDevPrior.cwiseProduct(trackOriginStdDevPrior).asDiagonal();

  // Alignment mask
  ActsAlignment::AlignmentMask alignmentMask =
      (ActsAlignment::AlignmentMask::Center1 |
       ActsAlignment::AlignmentMask::Center2 |
       ActsAlignment::AlignmentMask::Rotation2);

  // Alignment parameters solver

  LocalAlignmentParametersSolverConstraints::Config alignmentSolverCfg{};
  alignmentSolverCfg.alignmentMask = alignmentMask;
  LocalAlignmentParametersSolverConstraints alignmentSolver(alignmentSolverCfg,
                                                            logLevel);

  // LocalAlignmentParametersSolverSVD::Config alignmentSolverCfg{};
  // alignmentSolverCfg.alignmentMask = alignmentMask;
  // alignmentSolverCfg.maxSingularValueTol = 1e-4;
  // alignmentSolverCfg.singularValueGapTol = 9e-1;
  // LocalAlignmentParametersSolverSVD alignmentSolver(alignmentSolverCfg,
  //                                                   logLevel);

  // Alignment transform updater
  LocalAlignmentTransformUpdater::Config alignmentUpdaterCfg{};
  alignmentUpdaterCfg.alignmentStore = aStore;
  LocalAlignmentTransformUpdater alignmentUpdater(alignmentUpdaterCfg,
                                                  logLevel);

  // Number of refitting iterations
  std::size_t nRefittingIt = 1;

  // Annealing scheduler
  LinearAnnealingScheduler::Config annealingSchedulerCfg{};
  annealingSchedulerCfg.alphaStart = 1e3;
  annealingSchedulerCfg.alphaEnd = 1e3;
  annealingSchedulerCfg.nIt = nRefittingIt;

  auto annealingScheduler =
      std::make_shared<LinearAnnealingScheduler>(annealingSchedulerCfg);

  // GX2 fitter setup
  StraightLineGX2Fitter::Config gx2FitterCfg{};
  gx2FitterCfg.primaryIdx = goInst.primaryIdx;
  gx2FitterCfg.longIdx = goInst.longIdx;
  gx2FitterCfg.shortIdx = goInst.shortIdx;
  gx2FitterCfg.firstLayerGeoId = goInst.tcParameters.front().geoId;
  gx2FitterCfg.lastLayerGeoId = goInst.tcParameters.back().geoId;
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
      E320::E320TrackParametersEstimator::PropagationDirection::forward;

  auto trackParametersEstimator =
      std::make_shared<E320::E320TrackParametersEstimator>(
          trackParametersEstimatorCfg);

  // Alignment algorithm
  AlignmentAlgorithm::Config alignmentCfg{
      .inputSourceLinks = "SourceLinks",
      .inputTrackCandidates = "SeedsGuess",
      .inputTrackParameters = "TrackParametersGuess",
      .outputAlignmentParameters = "AlignmentParameters",
      .outputTrackParameters = "UpdatedTrackParameters",
      .alignmentFunction =
          AlignmentAlgorithm::makeAlignmentFunction(detector, field),
      .maxKFSteps = static_cast<std::size_t>(1e5),
      .kfExtensions = alignmentExtensions,
      .kfReferenceSurface = trackingRefSurface.get(),
      .chi2ONdfCutOff = 1e-16,
      .deltaChi2ONdfCutOff = {10, 1e-5},
      .maxAlignmentFitNumIt = 200,
      .alignmentMask = alignmentMask,
      .originCov = trackOriginCov,
      .constraints = {},
      .nRefittingIt = nRefittingIt,
      .annealingScheduler = annealingScheduler,
      .trackParametersEstimator = trackParametersEstimator};

  alignmentCfg.alignmentParametersSolver.connect<
      &LocalAlignmentParametersSolverConstraints::calculateAlignmentParameters>(
      &alignmentSolver);

  alignmentCfg.alignmentTransformUpdater
      .connect<&LocalAlignmentTransformUpdater::updateAlignmentParameters>(
          &alignmentUpdater);

  for (auto& det : detector->detectorElements()) {
    const auto& surface = det->surface();
    const auto& geoId = surface.geometryId().sensitive();
    if (surface.geometryId().sensitive() > goInst.tcParameters.front().geoId &&
        surface.geometryId().sensitive() <= goInst.tcParameters.back().geoId) {
      alignmentCfg.alignedDetElements.push_back(det.get());
    }
  }

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

  // Seed writer
  RootSeedWriter::Config seedWriterCfg;
  seedWriterCfg.inputSeeds = "SeedsGuess";
  seedWriterCfg.inputTrackParameters = "UpdatedTrackParameters";
  seedWriterCfg.inputSourceLinks = "SourceLinks";
  seedWriterCfg.treeName = "seeds";
  seedWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/data/"
      "seeds.root";

  sequencer.addWriter(
      std::make_shared<RootSeedWriter>(seedWriterCfg, logLevel));

  // Fitted track writer
  RootTrackWriter::Config trackWriterCfg;
  trackWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);
  trackWriterCfg.referenceSurface = trackingRefSurface.get();
  trackWriterCfg.inputTrackContainer = "TrackContainer";
  trackWriterCfg.inputTracks = "Tracks";
  trackWriterCfg.inputTrackParametersGuesses = "UpdatedTrackParameters";
  trackWriterCfg.treeName = "fitted-tracks";
  trackWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/data/"
      "fitted-tracks.root";

  sequencer.addWriter(
      std::make_shared<RootTrackWriter>(trackWriterCfg, logLevel));

  // Alignment parameters writer
  AlignmentParametersWriter::Config alignmentWriterCfg;
  alignmentWriterCfg.treeName = "alignment-parameters";
  alignmentWriterCfg.inputAlignmentResults = "AlignmentParameters";
  alignmentWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/data/"
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
