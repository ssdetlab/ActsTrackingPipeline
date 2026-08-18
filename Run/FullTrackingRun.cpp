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
#include "TrackingPipeline/EventData/MixedSourceLinkCalibrator.hpp"
#include "TrackingPipeline/EventData/MixedSourceLinkSurfaceAccessor.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryOptions.hpp"
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
      "FullTrackingRunConfig.toml";
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
  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> surfaceMap;
  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
      gx2FitterSurfaceMap;
  for (const auto& vol : detector->volumes()) {
    for (const auto& surf : vol->surfaces()) {
      const auto& geoId = surf->geometryId();
      if (geoId.sensitive() != 0u) {
        surfaceMap[geoId] = surf;
        if (geoId.sensitive() >= goInst.tcParameters.front().geoId &&
            geoId.sensitive() <= goInst.tcParameters.back().geoId) {
          gx2FitterSurfaceMap[geoId] = surf;
        }
      }
    }
  }

  // Initialize alignment store
  auto aStore = detail::makeAlignmentStore(gctx, detector.get());

  // AlignmentParametersProvider::Config alignmentProviderCfg;
  // alignmentProviderCfg.filePath =
  //     "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/data/"
  //     "alignment/local_2026_data/big_shifts/aligned/alignment-parameters.root";
  // alignmentProviderCfg.treeName = "alignment-parameters";
  // AlignmentParametersProvider alignmentProvider(alignmentProviderCfg);
  // aStore = alignmentProvider.getAlignmentStore();

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
  gctx = Acts::GeometryContext{alignCtx};

  // --------------------------------------------------------------
  // The magnetic field setup

  auto field = E320::buildMagField(gctx);

  E320::E320MagneticFieldParametersProvider::Config magFieldProviderCfg;
  magFieldProviderCfg.treeName =
      getEntryStr("E320MagneticFieldParametersProvider", "treeName");
  std::string inDirMagnets =
      getEntryStr("E320MagneticFieldParametersProvider", "inDirMagnets");

  // Get the paths to the files in the directory
  for (const auto& entry : std::filesystem::directory_iterator(inDirMagnets)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".root") {
      continue;
    }
    std::string pathToFile = entry.path();
    magFieldProviderCfg.filePaths.push_back(pathToFile);
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

  // Mixed surface accessor
  SimpleSourceLink::SurfaceAccessor simpleSurfaceAccessor{detector.get()};
  ExtendedSourceLink::SurfaceAccessor extendedSurfaceAccessor{detector.get()};
  MixedSourceLinkSurfaceAccessor mixedSurfaceAccessor;
  mixedSurfaceAccessor.connect<&SimpleSourceLink::SurfaceAccessor::operator(),
                               SimpleSourceLink>(&simpleSurfaceAccessor);
  mixedSurfaceAccessor.connect<&ExtendedSourceLink::SurfaceAccessor::operator(),
                               ExtendedSourceLink>(&extendedSurfaceAccessor);

  // Mixed calibrator
  MixedSourceLinkCalibrator<KFFitterTrajectory> mixedCalibrator;
  mixedCalibrator.connect<&simpleSourceLinkCalibrator<KFFitterTrajectory>,
                          SimpleSourceLink>();
  mixedCalibrator.connect<
      &extendedSourceLinkBackwardsPhiCorrectionCalibrator<KFFitterTrajectory>,
      ExtendedSourceLink>();

  // Setup the sequencer
  Sequencer::Config seqCfg;
  seqCfg.skip = getEntrySizeT("Sequencer", "skip");
  seqCfg.events = getEntrySizeT("Sequencer", "events");
  seqCfg.numThreads = getEntrySizeT("Sequencer", "numThreads");
  seqCfg.trackFpes = getEntryBool("Sequencer", "trackFpes");
  seqCfg.logLevel = logLevel;
  Sequencer sequencer(seqCfg);

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));
  sequencer.addContextDecorator(
      std::make_shared<MagneticFieldContextDecorator>(magFieldCollection));

  // Add the sim data reader
  E320::E320RootDataReader::Config readerCfg;
  readerCfg.outputSourceLinks =
      getEntryStr("E320RootDataReader", "outputSourceLinks");
  readerCfg.outputDetSourceLinkIndices =
      getEntryStr("E320RootDataReader", "outputDetSourceLinkIndices");
  readerCfg.outputBpmSourceLinkIndices =
      getEntryStr("E320RootDataReader", "outputBpmSourceLinkIndices");
  readerCfg.outputEventMetaData =
      getEntryStr("E320RootDataReader", "outputEventMetaData");
  readerCfg.treeName = getEntryStr("E320RootDataReader", "treeName");
  readerCfg.eventKey = getEntryStr("E320RootDataReader", "eventKey");
  readerCfg.requireEpicsParity =
      getEntryBool("E320RootDataReader", "requireEpicsParity");
  readerCfg.requiredEpicsParity = E320::E320RootDataReader::EpicsParity(
      getEntrySizeT("E320RootDataReader", "requiredEpicsParity"));
  readerCfg.maxOccupancy = getEntrySizeT("E320RootDataReader", "maxOccupancy");
  readerCfg.minGeoId = goInst.tcParameters.front().geoId;
  readerCfg.maxGeoId = goInst.tcParameters.back().geoId;
  readerCfg.surfaceMap = surfaceMap;

  std::string inDirClusters =
      getEntryStr("E320RootDataReader", "inDirClusters");

  // Get the paths to the files in the directory
  for (const auto& entry : std::filesystem::directory_iterator(inDirClusters)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".root") {
      continue;
    }
    std::string pathToFile = entry.path();
    readerCfg.filePaths.push_back(pathToFile);
  }

  // Add the reader to the sequencer
  sequencer.addReader(
      std::make_shared<E320::E320RootDataReader>(readerCfg, logLevel));

  // --------------------------------------------------------------
  // Seeding setup

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

  // Covariance prior
  Acts::BoundVector trackOriginStdDevPrior;
  trackOriginStdDevPrior[Acts::eBoundLoc0] = 10_mm;
  trackOriginStdDevPrior[Acts::eBoundLoc1] = 10_mm;
  trackOriginStdDevPrior[Acts::eBoundPhi] = 1_rad;
  trackOriginStdDevPrior[Acts::eBoundTheta] = 1_rad;
  trackOriginStdDevPrior[Acts::eBoundQOverP] = 1_e / 1_GeV;
  trackOriginStdDevPrior[Acts::eBoundTime] = 1_fs;

  // Track parameters estimator
  E320::E320TrackParametersEstimator::Config trackParametersEstimatorCfg{};
  trackParametersEstimatorCfg.gx2Fitter = gx2Fitter;
  trackParametersEstimatorCfg.referenceSurface = seedingRefSurface.get();
  trackParametersEstimatorCfg.originCov =
      trackOriginStdDevPrior.cwiseProduct(trackOriginStdDevPrior).asDiagonal();

  trackParametersEstimatorCfg.nIterations =
      getEntrySizeT("E320TrackParametersEstimator", "nIterations");
  trackParametersEstimatorCfg.maxChi2 =
      getEntryDouble("E320TrackParametersEstimator", "maxChi2");
  trackParametersEstimatorCfg.propDirection =
      E320::E320TrackParametersEstimator::PropagationDirection(
          getEntryInt("E320TrackParametersEstimator", "propDirection"));

  auto trackParametersEstimator =
      std::make_shared<E320::E320TrackParametersEstimator>(
          trackParametersEstimatorCfg);

  // HT seeder setup
  HoughTransformSeeder::Config htSeederCfg{};
  htSeederCfg.primaryIdx = goInst.primaryIdx;
  htSeederCfg.longIdx = goInst.longIdx;
  htSeederCfg.shortIdx = goInst.shortIdx;

  HoughTransformSeeder::Options htSeederOpt;
  htSeederOpt.nCellsThetaShort =
      getEntrySizeT("HoughTransformSeeder", "nCellsThetaShort");
  htSeederOpt.nCellsRhoShort =
      getEntrySizeT("HoughTransformSeeder", "nCellsRhoShort");
  htSeederOpt.nCellsThetaLong =
      getEntrySizeT("HoughTransformSeeder", "nCellsThetaLong");
  htSeederOpt.nCellsRhoLong =
      getEntrySizeT("HoughTransformSeeder", "nCellsRhoLong");
  htSeederOpt.minXCount = getEntrySizeT("HoughTransformSeeder", "minXCount");

  htSeederOpt.boundBoxCenterPrimary = goInst.tcCenterPrimary;
  htSeederOpt.boundBoxCenterLong = goInst.tcCenterLong;
  htSeederOpt.boundBoxCenterShort = goInst.tcCenterShort;

  htSeederOpt.boundBoxHalfPrimary = goInst.tcHalfPrimary;
  htSeederOpt.boundBoxHalfLong = goInst.tcHalfLong;
  htSeederOpt.boundBoxHalfShort = goInst.tcHalfShort;

  htSeederOpt.surfaceMap = surfaceMap;

  // Seeding algorithm
  E320::E320SeedingAlgorithm::Config seedingAlgoCfg;
  seedingAlgoCfg.inputSourceLinks =
      getEntryStr("E320SeedingAlgorithm", "inputSourceLinks");
  seedingAlgoCfg.inputDetSourceLinkIndices =
      getEntryStr("E320SeedingAlgorithm", "inputDetSourceLinkIndices");
  seedingAlgoCfg.inputBpmSourceLinkIndices =
      getEntryStr("E320SeedingAlgorithm", "inputBpmSourceLinkIndices");
  seedingAlgoCfg.outputSeeds =
      getEntryStr("E320SeedingAlgorithm", "outputSeeds");
  seedingAlgoCfg.outputTrackParameters =
      getEntryStr("E320SeedingAlgorithm", "outputTrackParameters");
  seedingAlgoCfg.minLayers = getEntrySizeT("E320SeedingAlgorithm", "minLayers");
  seedingAlgoCfg.maxLayers = getEntrySizeT("E320SeedingAlgorithm", "maxLayers");
  seedingAlgoCfg.htSeeder = std::make_shared<HoughTransformSeeder>(htSeederCfg);
  seedingAlgoCfg.htOptions = htSeederOpt;
  seedingAlgoCfg.trackParametersEstimator = trackParametersEstimator;

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
      .connect<&MixedSourceLinkCalibrator<KFFitterTrajectory>::operator()>(
          &mixedCalibrator);
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
      .connect<&MixedSourceLinkSurfaceAccessor::operator()>(
          &mixedSurfaceAccessor);

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
      .inputTrackCandidates =
          getEntryStr("KFTrackFittingAlgorithm", "inputTrackCandidates"),
      .inputTrackParameters =
          getEntryStr("KFTrackFittingAlgorithm", "inputTrackParameters"),
      .inputSourceLinks =
          getEntryStr("KFTrackFittingAlgorithm", "inputSourceLinks"),
      .outputTrackContainer =
          getEntryStr("KFTrackFittingAlgorithm", "outputTrackContainer"),
      .outputTracks = getEntryStr("KFTrackFittingAlgorithm", "outputTracks"),
      .fitter = fitter,
      .maxSteps = getEntrySizeT("KFTrackFittingAlgorithm", "maxSteps"),
      .kfExtensions = kfExtensions,
      .referenceSurface = trackingRefSurface.get()};

  sequencer.addAlgorithm(
      std::make_shared<KFTrackFittingAlgorithm>(fitterCfg, logLevel));

  // --------------------------------------------------------------
  // Event write out

  // Cluster writer
  RootMeasurementWriter::Config measurementWriterCfg{};
  measurementWriterCfg.inputSourceLinks =
      getEntryStr("RootMeasurementWriter", "inputSourceLinks");
  measurementWriterCfg.treeName =
      getEntryStr("RootMeasurementWriter", "treeName");
  measurementWriterCfg.filePath =
      getEntryStr("RootMeasurementWriter", "filePath");

  sequencer.addWriter(
      std::make_shared<RootMeasurementWriter>(measurementWriterCfg, logLevel));

  // Seed writer
  RootSeedWriter::Config seedWriterCfg;
  seedWriterCfg.inputSeeds = getEntryStr("RootSeedWriter", "inputSeeds");
  seedWriterCfg.inputTrackParameters =
      getEntryStr("RootSeedWriter", "inputTrackParameters");
  seedWriterCfg.inputSourceLinks =
      getEntryStr("RootSeedWriter", "inputSourceLinks");
  seedWriterCfg.treeName = getEntryStr("RootSeedWriter", "treeName");
  seedWriterCfg.filePath = getEntryStr("RootSeedWriter", "filePath");

  sequencer.addWriter(
      std::make_shared<RootSeedWriter>(seedWriterCfg, logLevel));

  // Fitted track writer
  E320::E320RootTrackWriter::Config trackWriterCfg;
  trackWriterCfg.surfaceAccessor
      .connect<&MixedSourceLinkSurfaceAccessor::operator()>(
          &mixedSurfaceAccessor);
  trackWriterCfg.referenceSurface = trackingRefSurface.get();
  trackWriterCfg.inputTrackContainer =
      getEntryStr("E320RootTrackWriter", "inputTrackContainer");
  trackWriterCfg.inputTracks =
      getEntryStr("E320RootTrackWriter", "inputTracks");
  trackWriterCfg.inputTrackParametersGuesses =
      getEntryStr("E320RootTrackWriter", "inputTrackParametersGuesses");
  trackWriterCfg.inputEventMetaData =
      getEntryStr("E320RootTrackWriter", "inputEventMetaData");
  trackWriterCfg.treeName = getEntryStr("E320RootTrackWriter", "treeName");
  trackWriterCfg.filePath = getEntryStr("E320RootTrackWriter", "filePath");

  sequencer.addWriter(
      std::make_shared<E320::E320RootTrackWriter>(trackWriterCfg, logLevel));

  // Magnetic field writer
  E320::E320MagneticFieldWriter::Config magFieldWriterCfg;
  magFieldWriterCfg.treeName =
      getEntryStr("E320MagneticFieldWriter", "treeName");
  magFieldWriterCfg.filePath =
      getEntryStr("E320MagneticFieldWriter", "filePath");

  sequencer.addWriter(std::make_shared<E320::E320MagneticFieldWriter>(
      magFieldWriterCfg, logLevel));

  return sequencer.run();
}
