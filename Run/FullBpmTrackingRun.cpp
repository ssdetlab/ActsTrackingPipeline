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
#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/Geometry/GeometryContext.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

#include <unistd.h>

#include "TrackingPipeline/Alignment/detail/AlignmentStoreBuilders.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Infrastructure/TypeDefinitions.hpp"
#include "TrackingPipeline/Io/AlignmentParametersProvider.hpp"
#include "TrackingPipeline/Io/E320RootDataReader.hpp"
#include "TrackingPipeline/Io/RootMeasurementWriter.hpp"
#include "TrackingPipeline/Io/RootSeedWriter.hpp"
#include "TrackingPipeline/Io/RootTrackWriter.hpp"
#include "TrackingPipeline/MagneticField/ConstantMagField.hpp"
#include "TrackingPipeline/MagneticField/IdealQuadrupoleMagField.hpp"
#include "TrackingPipeline/MagneticField/MagneticFieldContextDecorator.hpp"
#include "TrackingPipeline/MagneticField/MagneticFieldParametersContext.hpp"
#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"
#include "TrackingPipeline/TrackFinding/E320SeedingAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/E320TrackParametersEstimator.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"
#include "TrackingPipeline/TrackFinding/StraightLineTrackParametersEstimator.hpp"
#include "TrackingPipeline/TrackFinding/TryAllSeedingAlgorithm.hpp"
#include "TrackingPipeline/TrackFitting/KFTrackFittingAlgorithm.hpp"
#include "TrackingPipeline/TrackFitting/StraightLineGX2Fitter.hpp"
#include "TrackingPipeline/Utilities/ThetaMcsRmsCalculator.hpp"

using namespace Acts::UnitLiterals;

std::unique_ptr<const E320::GeometryOptions> E320::GeometryOptions::m_instance =
    nullptr;

int main() {
  // Geometry constraints instance
  const auto& goInst = *E320::GeometryOptions::instance();

  // Set the log level
  Acts::Logging::Level logLevel = Acts::Logging::INFO;

  // Initialize contexts
  Acts::GeometryContext gctx;

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
        if (surf->geometryId().sensitive() >= goInst.bpm0Parameters.geoId &&
            surf->geometryId().sensitive() <= goInst.bpm3Parameters.geoId) {
          gx2FitterSurfaceMap[surf->geometryId()] = surf;
        }
      }
    }
  }

  // Initialize alignment store
  auto aStore = detail::makeAlignmentStore(gctx, detector.get());

  AlignmentParametersProvider::Config alignmentProviderCfg;
  alignmentProviderCfg.filePath =
      // "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/data/"
      // "alignment/local_feb_2025_data/local_solid_material_uniform_errors/"
      // "bkg_features_isolation/sig/aligned/"
      // "alignment-parameters.root";
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/data/"
      "alignment/global_feb_2025_data/global_mapped_material_uniform_errors/"
      "aligned/"
      "alignment-parameters.root";
  alignmentProviderCfg.treeName = "alignment-parameters";
  AlignmentParametersProvider alignmentProvider(alignmentProviderCfg);
  aStore = alignmentProvider.getAlignmentStore();

  // Initialize alignment context
  AlignmentContext alignCtx(aStore);

  // Print alignment parameters
  Acts::GeometryContext defaultGctx;
  gctx = Acts::GeometryContext{alignCtx};
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive() != 0u) {
        std::cout << "-----------------------------------\n";
        std::cout << "SURFACE " << s->geometryId() << "\n";
        std::cout << "CENTER " << s->center(gctx).transpose() << " -- "
                  << s->center(defaultGctx).transpose() << "\n";
        std::cout << "DELTA "
                  << (s->center(gctx) - s->center(defaultGctx)).transpose() *
                         1e3
                  << "\n";
        std::cout << "NORMAL "
                  << s->normal(gctx, s->center(gctx), Acts::Vector3::UnitY())
                         .transpose()
                  << " -- "
                  << s->normal(gctx, s->center(defaultGctx),
                               Acts::Vector3::UnitY())
                         .transpose()
                  << "\n";
        std::cout << "ROTATION \n"
                  << s->transform(gctx).rotation() << " -- \n"
                  << "\n"
                  << s->transform(defaultGctx).rotation() << "\n";

        std::cout << "EXTENT "
                  << s->polyhedronRepresentation(gctx, 1000).extent()
                  << "\n -- \n"
                  << s->polyhedronRepresentation(defaultGctx, 1000).extent()
                  << "\n";
      }
    }
  }

  // // --------------------------------------------------------------
  // // The magnetic field setup

  // auto mStore1 = std::make_shared<MagneticFieldStore>();
  // mStore1->store = {
  //     {goInst.quad1Id,
  //      Acts::MagneticFieldProvider::Cache(
  //          std::in_place_type<IdealQuadrupoleMagField::Cache>, 0)},
  //     {goInst.quad2Id,
  //      Acts::MagneticFieldProvider::Cache(
  //          std::in_place_type<IdealQuadrupoleMagField::Cache>, 0)},
  //     {goInst.quad3Id,
  //      Acts::MagneticFieldProvider::Cache(
  //          std::in_place_type<IdealQuadrupoleMagField::Cache>, 0)},
  //     {goInst.xCorrectorId, Acts::MagneticFieldProvider::Cache(
  //                               std::in_place_type<ConstantMagField::Cache>,
  //                               Acts::Vector3(0, 0.026107_T, 0))},
  //     {goInst.dipoleId, Acts::MagneticFieldProvider::Cache(
  //                           std::in_place_type<ConstantMagField::Cache>,
  //                           Acts::Vector3(0, 0, -0.2192_T))}};

  // auto mStore2 = std::make_shared<MagneticFieldStore>();
  // mStore2->store = {
  //     {goInst.quad1Id,
  //      Acts::MagneticFieldProvider::Cache(
  //          std::in_place_type<IdealQuadrupoleMagField::Cache>, 0)},
  //     {goInst.quad2Id,
  //      Acts::MagneticFieldProvider::Cache(
  //          std::in_place_type<IdealQuadrupoleMagField::Cache>, 0)},
  //     {goInst.quad3Id,
  //      Acts::MagneticFieldProvider::Cache(
  //          std::in_place_type<IdealQuadrupoleMagField::Cache>, 0)},
  //     {goInst.xCorrectorId, Acts::MagneticFieldProvider::Cache(
  //                               std::in_place_type<ConstantMagField::Cache>,
  //                               Acts::Vector3(0, 0.026107_T, 0))},
  //     {goInst.dipoleId, Acts::MagneticFieldProvider::Cache(
  //                           std::in_place_type<ConstantMagField::Cache>,
  //                           Acts::Vector3(0, 0, -0.2192_T))}};

  // MagneticFieldStoreCollection mFieldStoreCollection = {{0, mStore1},
  //                                                       {1000, mStore2}};
  // auto mFieldParametersContext =
  //     std::make_shared<MagneticFieldParametersContext>(mFieldStoreCollection);

  auto field = E320::buildMagField(gctx);

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

  // Seeding reference surface
  Acts::Transform3 seedingRefSurfTransform = Acts::Transform3::Identity();
  seedingRefSurfTransform.translation() =
      Acts::Vector3(goInst.bpm0CenterPrimary - 0.1_mm, 0, 0);
  seedingRefSurfTransform.rotate(refSurfToWorldRotationX);
  seedingRefSurfTransform.rotate(refSurfToWorldRotationY);
  seedingRefSurfTransform.rotate(refSurfToWorldRotationZ);

  auto seedingRefSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      seedingRefSurfTransform,
      std::make_shared<Acts::RectangleBounds>(halfX, halfY));

  Acts::GeometryIdentifier seedingRefSurfaceGeoId;
  seedingRefSurfaceGeoId.setExtra(1);
  seedingRefSurface->assignGeometryId(seedingRefSurfaceGeoId);

  // Tracking reference surface
  Acts::Transform3 trackingRefSurfTransform = Acts::Transform3::Identity();
  trackingRefSurfTransform.translation() =
      Acts::Vector3(goInst.bpm0CenterPrimary - 0.1_mm, 0, 0);
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
  seqCfg.events = 100;
  seqCfg.numThreads = 32;
  seqCfg.trackFpes = false;
  seqCfg.logLevel = logLevel;
  Sequencer sequencer(seqCfg);

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));

  // Add the sim data reader
  E320::E320RootDataReader::Config readerCfg;
  readerCfg.outputSourceLinks = "SourceLinks";
  readerCfg.outputDetSourceLinkIndices = "DetSourceLinkIndices";
  readerCfg.outputBpmSourceLinkIndices = "BpmSourceLinkIndices";
  readerCfg.treeName = "MyTree";
  readerCfg.eventKey = "event";
  readerCfg.minGeoId = goInst.bpm0Parameters.geoId;
  readerCfg.maxGeoId = goInst.bpm3Parameters.geoId;
  readerCfg.surfaceMap = surfaceMap;
  std::string pathToDir =
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/E320Shift_Feb_2025/processed/"
      "data_Run502";

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
      std::make_shared<E320::E320RootDataReader>(readerCfg, logLevel));

  // --------------------------------------------------------------
  // Seeding setup

  // GX2 fitter setup
  StraightLineGX2Fitter::Config gx2FitterCfg{};
  gx2FitterCfg.primaryIdx = goInst.primaryIdx;
  gx2FitterCfg.longIdx = goInst.longIdx;
  gx2FitterCfg.shortIdx = goInst.shortIdx;
  gx2FitterCfg.firstLayerGeoId = goInst.bpm0Parameters.geoId;
  gx2FitterCfg.lastLayerGeoId = goInst.bpm3Parameters.geoId;
  gx2FitterCfg.surfaceMap = gx2FitterSurfaceMap;

  auto gx2Fitter = std::make_shared<StraightLineGX2Fitter>(gx2FitterCfg);

  // Covariance prior
  Acts::BoundVector trackOriginStdDevPrior;
  trackOriginStdDevPrior[Acts::eBoundLoc0] = 100_mm;
  trackOriginStdDevPrior[Acts::eBoundLoc1] = 100_mm;
  trackOriginStdDevPrior[Acts::eBoundPhi] = 10_degree;
  trackOriginStdDevPrior[Acts::eBoundTheta] = 10_degree;
  trackOriginStdDevPrior[Acts::eBoundQOverP] = 1 / 0.01_GeV;
  trackOriginStdDevPrior[Acts::eBoundTime] = 1_fs;

  // Track parameters estimator
  StraightLineTrackParametersEstimator::Config trackParametersEstimatorCfg{
      .gx2Fitter = gx2Fitter,
      .nIterations = 2,
      .maxChi2 = std::numeric_limits<double>::max(),
      .referenceSurface = seedingRefSurface.get(),
      .absMomentum = 10_GeV,
      .particleHypothesis = Acts::ParticleHypothesis::electron(),
      .charge = -1_e,
      .originCov = trackOriginStdDevPrior.cwiseProduct(trackOriginStdDevPrior)
                       .asDiagonal(),
      .thetaRms = 0};

  auto trackParametersEstimator =
      std::make_shared<StraightLineTrackParametersEstimator>(
          trackParametersEstimatorCfg);

  // Seeding algorithm
  TryAllSeedingAlgorithm::Config seedingAlgoCfg;
  seedingAlgoCfg.inputSourceLinks = "SourceLinks";
  seedingAlgoCfg.inputSourceLinkIndices = "BpmSourceLinkIndices";
  seedingAlgoCfg.outputSeeds = "Seeds";
  seedingAlgoCfg.outputTrackParameters = "TrackParameters";
  seedingAlgoCfg.trackParametersEstimator = trackParametersEstimator;
  seedingAlgoCfg.minLayers = 4;
  seedingAlgoCfg.maxLayers = 4;

  sequencer.addAlgorithm(
      std::make_shared<TryAllSeedingAlgorithm>(seedingAlgoCfg, logLevel));

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
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/data/"
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
  trackWriterCfg.inputTrackParametersGuesses = "TrackParameters";
  trackWriterCfg.treeName = "fitted-tracks";
  trackWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/data/"
      "fitted-tracks.root";

  sequencer.addWriter(
      std::make_shared<RootTrackWriter>(trackWriterCfg, logLevel));

  return sequencer.run();
}
