#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>

#include <unistd.h>

#include "TrackingPipeline/Geometry/ApollonGeometry.hpp"
#include "TrackingPipeline/Geometry/ApollonGeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/AlignmentParametersProvider.hpp"
#include "TrackingPipeline/Io/ApollonRootSimDataReader.hpp"
#include "TrackingPipeline/Io/RootSimClusterReader.hpp"
#include "TrackingPipeline/Io/RootSimSeedWriter.hpp"
#include "TrackingPipeline/Io/RootSimTrackWriter.hpp"
#include "TrackingPipeline/TrackFinding/ApollonSeedingAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"
#include "TrackingPipeline/TrackFitting/KFTrackFittingAlgorithm.hpp"

// Propagator short-hands
using ActionList = Acts::ActionList<>;
using AbortList = Acts::AbortList<Acts::EndOfWorldReached>;

using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                    Acts::Experimental::DetectorNavigator>;
using PropagatorOptions =
    typename Propagator::template Options<ActionList, AbortList>;

// KF short-hands
using RecoTrajectory = KFTrackFittingAlgorithm::Trajectory;
using RecoTrackContainer = KFTrackFittingAlgorithm::TrackContainer;
using KF = Acts::KalmanFitter<Propagator, RecoTrajectory>;

using namespace Acts::UnitLiterals;

namespace ag = ApollonGeometry;

std::unique_ptr<const ag::GeometryOptions> ag::GeometryOptions::m_instance =
    nullptr;

int main() {
  const auto& goInst = *ag::GeometryOptions::instance();

  // Set the log level
  Acts::Logging::Level logLevel = Acts::Logging::DEBUG;

  // Dummy context and options
  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext cctx;

  // --------------------------------------------------------------
  // Detector setup

  auto detector = ApollonGeometry::buildDetector(gctx);

  std::map<Acts::GeometryIdentifier, const Acts::Surface*> surfaceMap;
  for (const auto& vol : detector->volumes()) {
    std::cout << "------------------------------------------\n";
    std::cout << vol->name() << "\n";
    std::cout << vol->extent(gctx);
    std::cout << "Surfaces:\n";
    for (const auto& surf : vol->surfaces()) {
      std::cout << surf->geometryId() << "\n";
      std::cout << surf->center(gctx) << "\n";
      std::cout << surf->polyhedronRepresentation(gctx, 1000).extent() << "\n";
      if (surf->geometryId().sensitive()) {
        surfaceMap[surf->geometryId()] = surf;
      }
    }
  }

  AlignmentParametersProvider::Config alignmentProviderCfg1;
  alignmentProviderCfg1.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/fast_sim_data/"
      "alignment/det1/aligned/"
      "alignment-parameters.root";
  alignmentProviderCfg1.treeName = "alignment-parameters";
  AlignmentParametersProvider alignmentProvider1(alignmentProviderCfg1);
  auto aStore1 = alignmentProvider1.getAlignmentStore();

  AlignmentParametersProvider::Config alignmentProviderCfg2;
  alignmentProviderCfg2.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/fast_sim_data/"
      "alignment/det2/aligned/"
      "alignment-parameters.root";
  alignmentProviderCfg2.treeName = "alignment-parameters";
  AlignmentParametersProvider alignmentProvider2(alignmentProviderCfg2);
  auto aStore2 = alignmentProvider2.getAlignmentStore();

  auto aStore = std::make_shared<AlignmentContext::AlignmentStore>();
  for (const auto& entry : *aStore1) {
    aStore->insert(entry);
  }
  for (const auto& entry : *aStore2) {
    aStore->insert(entry);
  }

  // AlignmentParametersProvider::Config alignmentProviderCfg;
  // alignmentProviderCfg.filePath =
  //     "/home/romanurmanov/work/Apollon/tracking/out_data/fast_sim_data/"
  //     "alignment/global/aligned/"
  //     "alignment-parameters.root";
  // alignmentProviderCfg.treeName = "alignment-parameters";
  // AlignmentParametersProvider alignmentProvider(alignmentProviderCfg);
  // auto aStore = alignmentProvider.getAlignmentStore();

  AlignmentContext alignCtx(aStore);
  Acts::GeometryContext testCtx{alignCtx};
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive() && s->geometryId().sensitive() >= 20) {
        std::cout << "-----------------------------------\n";
        std::cout << "SURFACE " << s->geometryId() << "\n";
        std::cout << "CENTER " << s->center(testCtx).transpose() << " -- "
                  << s->center(Acts::GeometryContext()).transpose() << "\n";
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

  auto field = ApollonGeometry::buildMagField(gctx);

  // --------------------------------------------------------------
  // Event reading
  SimpleSourceLink::SurfaceAccessor surfaceAccessor{detector.get()};

  // Setup the sequencer
  Sequencer::Config seqCfg;
  // seqCfg.events = 1e1;
  // seqCfg.skip = 1;
  seqCfg.numThreads = 1;
  seqCfg.trackFpes = false;
  seqCfg.logLevel = logLevel;
  Sequencer sequencer(seqCfg);

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));

  // Add the sim data reader
  RootSimClusterReader::Config readerCfg;
  readerCfg.outputSourceLinks = "Measurements";
  readerCfg.outputSimClusters = "SimClusters";
  readerCfg.treeName = "clusters";
  readerCfg.minGeoId = 10;
  readerCfg.maxGeoId = 28;
  // std::string pathToDir =
  //     "/home/romanurmanov/work/Apollon/tracking/out_data/fast_sim_data/"
  //     "alignment/data";
  std::string pathToDir =
      "/home/romanurmanov/work/Apollon/tracking/out_data/fast_sim_data/"
      "alignment/data_test";

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
      std::make_shared<RootSimClusterReader>(readerCfg, logLevel));

  // // Add the sim data reader
  // ApollonIo::ApollonRootSimDataReader::Config readerCfg;
  // readerCfg.outputSourceLinks = "Measurements";
  // readerCfg.outputSimClusters = "SimClusters";
  // readerCfg.treeName = "clusters";
  // readerCfg.surfaceMap = surfaceMap;
  // readerCfg.eventSplit = false;
  // std::string pathToDir =
  //     "/home/romanurmanov/work/Apollon/geant4_sims/g4_clustering/out_data/"
  //     "test_clusters";

  // // Get the paths to the files in the directory
  // for (const auto& entry : std::filesystem::directory_iterator(pathToDir)) {
  //   if (!entry.is_regular_file() || entry.path().extension() != ".root") {
  //     continue;
  //   }
  //   std::string pathToFile = entry.path();
  //   readerCfg.filePaths.push_back(pathToFile);
  // }

  // // Add the reader to the sequencer
  // sequencer.addReader(std::make_shared<ApollonIo::ApollonRootSimDataReader>(
  //     readerCfg, logLevel));

  // --------------------------------------------------------------
  // HT seeding setup

  HoughTransformSeeder::Config htSeederCfg;
  htSeederCfg.boundBoxHalfX =
      goInst.tc1HalfPrimary - goInst.chipVolumeHalfSpacing * 2;
  htSeederCfg.boundBoxHalfY = goInst.tcHalfLong;
  htSeederCfg.boundBoxHalfZ = goInst.tcHalfShort;

  htSeederCfg.nCellsThetaX = 500;
  htSeederCfg.nCellsRhoX = 4000;

  htSeederCfg.nCellsThetaY = 500;
  htSeederCfg.nCellsRhoY = 4000;

  htSeederCfg.minSeedSize = 5;
  htSeederCfg.maxSeedSize = 1000;

  htSeederCfg.nLSIterations = 2;

  ApollonSeedingAlgorithm::Config seedingAlgoCfg;
  seedingAlgoCfg.htSeeder = std::make_shared<HoughTransformSeeder>(htSeederCfg);
  seedingAlgoCfg.inputSourceLinks = "Measurements";
  seedingAlgoCfg.outputSeeds = "Seeds";
  seedingAlgoCfg.minLayers = 10;
  seedingAlgoCfg.maxLayers = 10;
  seedingAlgoCfg.minSeedSize = 10;
  seedingAlgoCfg.maxSeedSize = 100;
  seedingAlgoCfg.minScanEnergy = 0.7_GeV;
  seedingAlgoCfg.maxScanEnergy = 0.7_GeV;
  seedingAlgoCfg.energyScanStep = 0.001_GeV;
  seedingAlgoCfg.maxConnectionDistance = 1.0;
  seedingAlgoCfg.scope = ApollonSeedingAlgorithm::SeedingScope::fullDetector;

  sequencer.addAlgorithm(
      std::make_shared<ApollonSeedingAlgorithm>(seedingAlgoCfg, logLevel));

  // --------------------------------------------------------------
  // Track fitting

  Acts::GainMatrixUpdater kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;

  // Initialize track fitter options
  Acts::KalmanFitterExtensions<RecoTrajectory> extensions;
  // Add calibrator
  extensions.calibrator.connect<&simpleSourceLinkCalibrator<RecoTrajectory>>();
  // Add the updater
  extensions.updater
      .connect<&Acts::GainMatrixUpdater::operator()<RecoTrajectory>>(
          &kfUpdater);
  // Add the smoother
  extensions.smoother
      .connect<&Acts::GainMatrixSmoother::operator()<RecoTrajectory>>(
          &kfSmoother);
  // Add the surface accessor
  extensions.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  auto propOptions = PropagatorOptions(gctx, mctx);

  propOptions.maxSteps = 1e5;

  auto options =
      Acts::KalmanFitterOptions(gctx, mctx, cctx, extensions, propOptions);

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

  Acts::Transform3 transform = Acts::Transform3::Identity();
  transform.rotate(refSurfToWorldRotationX);
  transform.rotate(refSurfToWorldRotationY);
  transform.rotate(refSurfToWorldRotationZ);

  auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      transform, std::make_shared<Acts::RectangleBounds>(halfX, halfY));

  Acts::GeometryIdentifier geoId;
  geoId.setExtra(1);
  refSurface->assignGeometryId(std::move(geoId));

  options.referenceSurface = refSurface.get();

  Acts::Experimental::DetectorNavigator::Config cfg;
  cfg.detector = detector.get();
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Experimental::DetectorNavigator kfNavigator(
      cfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));

  Acts::EigenStepper<> kfStepper(std::move(field));
  auto kfPropagator =
      Propagator(std::move(kfStepper), std::move(kfNavigator),
                 Acts::getDefaultLogger("Propagator", logLevel));

  const auto fitter = KF(
      kfPropagator, Acts::getDefaultLogger("DetectorKalmanFilter", logLevel));

  // Add the track fitting algorithm to the sequencer
  KFTrackFittingAlgorithm::Config fitterCfg{.inputTrackCandidates = "Seeds",
                                            .outputTracks = "Tracks",
                                            .fitter = fitter,
                                            .kfOptions = options};

  sequencer.addAlgorithm(
      std::make_shared<KFTrackFittingAlgorithm>(fitterCfg, logLevel));

  // --------------------------------------------------------------
  // Event write out

  // Seed writer
  auto seedWriterCfg = RootSimSeedWriter::Config();
  seedWriterCfg.inputSeeds = "Seeds";
  seedWriterCfg.inputTruthClusters = "SimClusters";
  seedWriterCfg.treeName = "seeds";
  seedWriterCfg.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/fast_sim_data/"
      "alignment/"
      "seeds.root";

  sequencer.addWriter(
      std::make_shared<RootSimSeedWriter>(seedWriterCfg, logLevel));

  // Fitted track writer
  auto trackWriterCfg = RootSimTrackWriter::Config();
  trackWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);
  trackWriterCfg.referenceSurface = refSurface.get();
  trackWriterCfg.inputTracks = "Tracks";
  trackWriterCfg.inputSimClusters = "SimClusters";
  trackWriterCfg.treeName = "fitted-tracks";
  trackWriterCfg.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/fast_sim_data/"
      "alignment/"
      "fitted-tracks.root";

  sequencer.addWriter(
      std::make_shared<RootSimTrackWriter>(trackWriterCfg, logLevel));

  return sequencer.run();
}
