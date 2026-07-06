#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
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

#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>

#include <unistd.h>

#include "TrackingPipeline/Alignment/detail/AlignmentStoreBuilders.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Infrastructure/TypeDefinitions.hpp"
#include "TrackingPipeline/Io/AlignmentParametersProvider.hpp"
#include "TrackingPipeline/Io/E320MagneticFieldParametersProvider.hpp"
#include "TrackingPipeline/Io/E320MagneticFieldWriter.hpp"
#include "TrackingPipeline/Io/E320RootDataReader.hpp"
#include "TrackingPipeline/Io/E320RootTrackWriter.hpp"
#include "TrackingPipeline/Io/RootMeasurementWriter.hpp"
#include "TrackingPipeline/Io/RootSeedWriter.hpp"
#include "TrackingPipeline/MagneticField/MagneticFieldContextDecorator.hpp"
#include "TrackingPipeline/TrackFinding/E320SeedingAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/E320TrackParametersEstimator.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"
#include "TrackingPipeline/TrackFitting/KFTrackFittingAlgorithm.hpp"
#include "TrackingPipeline/TrackFitting/StraightLineGX2Fitter.hpp"

using namespace Acts::UnitLiterals;

std::unique_ptr<const E320::GeometryOptions> E320::GeometryOptions::m_instance =
    nullptr;

int main() {
  const auto& goInst = *E320::GeometryOptions::instance();

  // Set the log level
  Acts::Logging::Level logLevel = Acts::Logging::INFO;

  // Dummy context and options
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

  auto detector = E320::buildDetector(gctx, materialDecorator);

  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> surfaceMap;
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
      if (surf->geometryId().sensitive() != 0u) {
        surfaceMap[surf->geometryId()] = surf;
        if (surf->geometryId().sensitive() >=
                goInst.tcParameters.front().geoId &&
            surf->geometryId().sensitive() <=
                goInst.tcParameters.back().geoId) {
          gx2FitterSurfaceMap[surf->geometryId()] = surf;
        }
      }
    }
  }

  auto aStore = detail::makeAlignmentStore(gctx, detector.get());

  AlignmentParametersProvider::Config alignmentProviderCfg;
  alignmentProviderCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/data/"
      "alignment/local_march_2026_data/sig/feb_data_procedure/test/aligned/"
      "alignment-parameters.root";
  alignmentProviderCfg.treeName = "alignment-parameters";
  AlignmentParametersProvider alignmentProvider(alignmentProviderCfg);
  aStore = alignmentProvider.getAlignmentStore();

  AlignmentContext alignCtx(aStore);
  Acts::GeometryContext testCtx{alignCtx};
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive() != 0u) {
        std::cout << "-----------------------------------\n";
        std::cout << "SURFACE " << s->geometryId() << "\n";
        std::cout << "CENTER " << s->center(testCtx).transpose() << " -- "
                  << s->center(Acts::GeometryContext()).transpose() << "\n";
        std::cout << "DELTA "
                  << (s->center(testCtx) - s->center(Acts::GeometryContext()))
                             .transpose() *
                         1e3
                  << "\n";

        std::cout << "NORMAL "
                  << s->normal(testCtx, s->center(testCtx),
                               Acts::Vector3::UnitY())
                         .transpose()
                  << " -- "
                  << s->normal(testCtx, s->center(Acts::GeometryContext()),
                               Acts::Vector3::UnitY())
                         .transpose()
                  << "\n";
        std::cout << "ROTATION \n"
                  << s->transform(testCtx).rotation() << " -- \n"
                  << "\n"
                  << s->transform(Acts::GeometryContext()).rotation() << "\n";

        std::cout << "EXTENT "
                  << s->polyhedronRepresentation(testCtx, 1000).extent()
                  << "\n -- \n"
                  << s->polyhedronRepresentation(Acts::GeometryContext(), 1000)
                         .extent()
                  << "\n";
      }
    }
  }
  std::cout << "ALIGNMENT COVARIANCE:\n" << aStore->covariance << "\n";
  gctx = Acts::GeometryContext{alignCtx};

  // --------------------------------------------------------------
  // The magnetic field setup

  auto field = E320::buildMagField(gctx);

  E320::E320MagneticFieldParametersProvider::Config magFieldProviderCfg;
  magFieldProviderCfg.treeName = "magnets";
  std::vector<std::string> inDirs = {
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/E320Shift_June_2026_positrons/processed/"
      "positron_runs/june_23/combined"};
  // "/home/romanurmanov/work/E320/E320Prototype/"
  // "E320Prototype_dataInRootFormat/E320Shift_June_2026_positrons/processed/"
  // "positron_runs/june_23/even"};

  // Get the paths to the files in the directory
  for (const auto& dir : inDirs) {
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
      if (!entry.is_regular_file() || entry.path().extension() != ".root") {
        continue;
      }
      std::string pathToFile = entry.path();
      magFieldProviderCfg.filePaths.push_back(pathToFile);
    }
  }
  E320::E320MagneticFieldParametersProvider magFieldProvider(
      magFieldProviderCfg);

  auto magFieldCollection = magFieldProvider.getMagneticFieldStoreCollection();

  // --------------------------------------------------------------
  // Reference surfaces
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

  Acts::Transform3 seedingRefSurfTransform = Acts::Transform3::Identity();
  seedingRefSurfTransform.translation() =
      Acts::Vector3(goInst.ipTcDistance - 0.1_mm, 0, 0);
  // Acts::Vector3(goInst.ipTcDistance + 2 * goInst.tcHalfPrimary + 0.1_mm, 0,
  //               0);
  seedingRefSurfTransform.rotate(refSurfToWorldRotationX);
  seedingRefSurfTransform.rotate(refSurfToWorldRotationY);
  seedingRefSurfTransform.rotate(refSurfToWorldRotationZ);

  auto seedingRefSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      seedingRefSurfTransform,
      std::make_shared<Acts::RectangleBounds>(halfX, halfY));

  Acts::GeometryIdentifier seedingRefSurfaceGeoId;
  seedingRefSurfaceGeoId.setExtra(1);
  seedingRefSurface->assignGeometryId(seedingRefSurfaceGeoId);

  Acts::Transform3 trackingRefSurfTransform = Acts::Transform3::Identity();
  trackingRefSurfTransform.translation() =
      // Acts::Vector3(goInst.beWindowCenterPrimary, 0, 0);
      // Acts::Vector3(goInst.bpm3CenterPrimary + 0.1_mm, 0, 0);
      Acts::Vector3(
          goInst.dipoleCenterPrimary + goInst.dipoleHalfPrimary + 0.01_mm, 0,
          0);
  trackingRefSurfTransform.rotate(refSurfToWorldRotationX);
  trackingRefSurfTransform.rotate(refSurfToWorldRotationY);
  trackingRefSurfTransform.rotate(refSurfToWorldRotationZ);

  auto trackingRefSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      trackingRefSurfTransform,
      std::make_shared<Acts::RectangleBounds>(halfX, halfY));

  Acts::GeometryIdentifier trackingRefSurfaceGeoId;
  trackingRefSurfaceGeoId.setExtra(2);
  trackingRefSurface->assignGeometryId(trackingRefSurfaceGeoId);

  // --------------------------------------------------------------
  // Event reading
  SimpleSourceLink::SurfaceAccessor surfaceAccessor{detector.get()};

  // Setup the sequencer
  Sequencer::Config seqCfg;
  // seqCfg.events = 10;
  // seqCfg.skip = 57500;
  seqCfg.numThreads = 1;
  seqCfg.trackFpes = false;
  seqCfg.logLevel = logLevel;
  Sequencer sequencer(seqCfg);

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));
  sequencer.addContextDecorator(
      std::make_shared<MagneticFieldContextDecorator>(magFieldCollection));

  // Add the sim data reader
  E320::E320RootDataReader::Config readerCfg;
  readerCfg.outputSourceLinks = "SourceLinks";
  readerCfg.outputDetSourceLinkIndices = "DetSourceLinkIndices";
  readerCfg.outputBpmSourceLinkIndices = "BpmSourceLinkIndices";
  readerCfg.outputEventMetaData = "EventMetaData";
  readerCfg.treeName = "data";
  readerCfg.eventKey = "event";
  readerCfg.minGeoId = goInst.tcParameters.front().geoId;
  readerCfg.maxGeoId = goInst.tcParameters.back().geoId;
  readerCfg.surfaceMap = surfaceMap;

  // Get the paths to the files in the directory
  for (const auto& dir : inDirs) {
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
      if (!entry.is_regular_file() || entry.path().extension() != ".root") {
        continue;
      }
      std::string pathToFile = entry.path();
      readerCfg.filePaths.push_back(pathToFile);
    }
  }

  // Add the reader to the sequencer
  sequencer.addReader(
      std::make_shared<E320::E320RootDataReader>(readerCfg, logLevel));

  // --------------------------------------------------------------
  // Seeding setup

  // GX2 fitter setup
  StraightLineGX2Fitter::Config gx2FitterCfg{};
  gx2FitterCfg.primaryIdx = goInst.primaryIdx;
  gx2FitterCfg.longIdx = goInst.longIdx;
  gx2FitterCfg.shortIdx = goInst.shortIdx;
  gx2FitterCfg.firstLayerGeoId = goInst.tcParameters.front().geoId;
  gx2FitterCfg.lastLayerGeoId = goInst.tcParameters.back().geoId;
  gx2FitterCfg.surfaceMap = gx2FitterSurfaceMap;

  auto gx2Fitter = std::make_shared<StraightLineGX2Fitter>(gx2FitterCfg);

  // HT seeder setup
  HoughTransformSeeder::Config htSeederCfg{};
  htSeederCfg.primaryIdx = goInst.primaryIdx;
  htSeederCfg.longIdx = goInst.longIdx;
  htSeederCfg.shortIdx = goInst.shortIdx;

  HoughTransformSeeder::Options htSeederOpt;
  htSeederOpt.boundBoxCenterPrimary = goInst.tcCenterPrimary;
  htSeederOpt.boundBoxCenterLong = goInst.tcCenterLong;
  htSeederOpt.boundBoxCenterShort = goInst.tcCenterShort;

  htSeederOpt.boundBoxHalfPrimary = goInst.tcHalfPrimary;
  htSeederOpt.boundBoxHalfLong = goInst.tcHalfLong;
  htSeederOpt.boundBoxHalfShort = goInst.tcHalfShort;

  htSeederOpt.nCellsThetaShort = 500;
  htSeederOpt.nCellsRhoShort = 2000;

  htSeederOpt.nCellsThetaLong = 500;
  htSeederOpt.nCellsRhoLong = 2000;

  htSeederOpt.surfaceMap = surfaceMap;

  htSeederOpt.minXCount = 5;

  // Covariance prior
  Acts::BoundVector trackOriginStdDevPrior;
  trackOriginStdDevPrior[Acts::eBoundLoc0] = 100_mm;
  trackOriginStdDevPrior[Acts::eBoundLoc1] = 100_mm;
  trackOriginStdDevPrior[Acts::eBoundPhi] = 10_degree;
  trackOriginStdDevPrior[Acts::eBoundTheta] = 10_degree;
  trackOriginStdDevPrior[Acts::eBoundQOverP] = 1 / 0.01_GeV;
  trackOriginStdDevPrior[Acts::eBoundTime] = 1_fs;

  // Track parameters estimator
  E320::E320TrackParametersEstimator::Config trackParametersEstimatorCfg{};
  trackParametersEstimatorCfg.gx2Fitter = gx2Fitter;
  trackParametersEstimatorCfg.nIterations = 2;
  trackParametersEstimatorCfg.maxChi2 = std::numeric_limits<double>::max();
  trackParametersEstimatorCfg.referenceSurface = seedingRefSurface.get();
  trackParametersEstimatorCfg.originCov =
      trackOriginStdDevPrior.cwiseProduct(trackOriginStdDevPrior).asDiagonal();
  trackParametersEstimatorCfg.propDirection =
      E320::E320TrackParametersEstimator::PropagationDirection::forward;

  auto trackParametersEstimator =
      std::make_shared<E320::E320TrackParametersEstimator>(
          trackParametersEstimatorCfg);

  // Seeding algorithm
  E320::E320SeedingAlgorithm::Config seedingAlgoCfg;
  seedingAlgoCfg.htSeeder = std::make_shared<HoughTransformSeeder>(htSeederCfg);
  seedingAlgoCfg.htOptions = htSeederOpt;
  seedingAlgoCfg.inputSourceLinks = "SourceLinks";
  seedingAlgoCfg.inputDetSourceLinkIndices = "DetSourceLinkIndices";
  seedingAlgoCfg.inputBpmSourceLinkIndices = "BpmSourceLinkIndices";
  seedingAlgoCfg.outputSeeds = "Seeds";
  seedingAlgoCfg.outputTrackParameters = "TrackParameters";
  seedingAlgoCfg.trackParametersEstimator = trackParametersEstimator;
  seedingAlgoCfg.minLayers = goInst.tcParameters.size();
  seedingAlgoCfg.maxLayers = goInst.tcParameters.size();

  sequencer.addAlgorithm(
      std::make_shared<E320::E320SeedingAlgorithm>(seedingAlgoCfg, logLevel));

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
      .inputTrackCandidates = "Seeds",
      .inputTrackParameters = "TrackParameters",
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

  // Cluster writer
  RootMeasurementWriter::Config measurementWriterCfg{};
  measurementWriterCfg.inputSourceLinks = "SourceLinks";
  measurementWriterCfg.treeName = "measurements";
  measurementWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/"
      "data/"
      "measurements.root";

  sequencer.addWriter(
      std::make_shared<RootMeasurementWriter>(measurementWriterCfg, logLevel));

  // Seed writer
  RootSeedWriter::Config seedWriterCfg;
  seedWriterCfg.inputSeeds = "Seeds";
  seedWriterCfg.inputTrackParameters = "TrackParameters";
  seedWriterCfg.inputSourceLinks = "SourceLinks";
  seedWriterCfg.treeName = "seeds";
  seedWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/"
      "data/"
      "seeds.root";

  sequencer.addWriter(
      std::make_shared<RootSeedWriter>(seedWriterCfg, logLevel));

  // Fitted track writer
  E320::E320RootTrackWriter::Config trackWriterCfg;
  trackWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);
  trackWriterCfg.referenceSurface = trackingRefSurface.get();
  trackWriterCfg.inputTrackContainer = "TrackContainer";
  trackWriterCfg.inputTracks = "Tracks";
  trackWriterCfg.inputTrackParametersGuesses = "TrackParameters";
  trackWriterCfg.inputEventMetaData = "EventMetaData";
  trackWriterCfg.treeName = "fitted-tracks";
  trackWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/"
      "data/"
      "fitted-tracks.root";

  sequencer.addWriter(
      std::make_shared<E320::E320RootTrackWriter>(trackWriterCfg, logLevel));

  // Magnetic field writer
  auto magFieldWriterCfg = E320::E320MagneticFieldWriter::Config();
  magFieldWriterCfg.treeName = "magnets";
  magFieldWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/"
      "data/"
      "magnets.root";

  sequencer.addWriter(std::make_shared<E320::E320MagneticFieldWriter>(
      magFieldWriterCfg, logLevel));

  return sequencer.run();
}
