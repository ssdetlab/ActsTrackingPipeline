#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>

#include <unistd.h>

#include "TrackingPipeline/Geometry/ApollonGeometry.hpp"
#include "TrackingPipeline/Geometry/ApollonGeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/ApollonRootSimDataReader.hpp"
#include "TrackingPipeline/Io/RootSimClusterReader.hpp"
#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"
#include "TrackingPipeline/Io/RootSimSeedWriter.hpp"
#include "TrackingPipeline/Io/RootSimTrackCandidateWriter.hpp"
#include "TrackingPipeline/Io/RootSimTrackWriter.hpp"
#include "TrackingPipeline/TrackFinding/ApollonSeedingAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/CKFTrackFindingAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"
#include "TrackingPipeline/TrackFinding/IdealSeedingAlgorithm.hpp"
#include "TrackingPipeline/TrackFitting/KFTrackFittingAlgorithm.hpp"

// Propagator short-hands
using ActionList = Acts::ActionList<>;
using AbortList = Acts::AbortList<Acts::EndOfWorldReached>;

using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                    Acts::Experimental::DetectorNavigator>;
using PropagatorOptions =
    typename Propagator::template Options<ActionList, AbortList>;

// CKF short-hands
using CKFTrackContainer = CKFTrackFindingAlgorithm::TrackContainer;

using TrackStateContainerBackend =
    CKFTrackFindingAlgorithm::TrackStateContainerBackend;

// KF short-hands
using RecoTrajectory = KFTrackFittingAlgorithm::Trajectory;
using RecoTrackContainer = KFTrackFittingAlgorithm::TrackContainer;
using KF = Acts::KalmanFitter<Propagator, RecoTrajectory>;

using namespace Acts::UnitLiterals;

namespace ag = ApollonGeometry;

std::unique_ptr<const ag::GeometryOptions> ag::GeometryOptions::m_instance =
    nullptr;

int main() {
  // Set the log level
  Acts::Logging::Level logLevel = Acts::Logging::DEBUG;

  // Dummy context and options
  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext cctx;

  // --------------------------------------------------------------
  // Detector setup

  auto detector = ApollonGeometry::buildDetector(gctx);

  for (const auto& vol : detector->volumes()) {
    std::cout << "------------------------------------------\n";
    std::cout << vol->name() << "\n";
    std::cout << vol->extent(gctx);
    std::cout << "Surfaces:\n";
    for (const auto& surf : vol->surfaces()) {
      std::cout << surf->geometryId() << "\n";
      std::cout << surf->polyhedronRepresentation(gctx, 1000).extent() << "\n";
    }
  }

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

  // Add the sim data reader
  ApollonIo::ApollonRootSimDataReader::Config readerCfg;
  readerCfg.outputSourceLinks = "Measurements";
  readerCfg.outputSimClusters = "SimClusters";
  readerCfg.treeName = "particles";
  readerCfg.eventSplit = false;
  std::string pathToDir =
      "/home/romanurmanov/work/Apollon/geant4_sims/al_window_flange/out_data/"
      "test";

  // Get the paths to the files in the directory
  for (const auto& entry : std::filesystem::directory_iterator(pathToDir)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".root") {
      continue;
    }
    std::string pathToFile = entry.path();
    readerCfg.filePaths.push_back(pathToFile);
  }

  // Add the reader to the sequencer
  sequencer.addReader(std::make_shared<ApollonIo::ApollonRootSimDataReader>(
      readerCfg, logLevel));

  // // Add the sim data reader
  // RootSimClusterReader::Config readerCfg;
  // readerCfg.outputSourceLinks = "Measurements";
  // readerCfg.outputSimClusters = "SimClusters";
  // readerCfg.treeName = "clusters";
  // std::string pathToDir =
  //     "/home/romanurmanov/work/Apollon/tracking/out_data/fast_sim_data/"
  //     "sig_plus_bkg";

  // // Get the paths to the files in the directory
  // for (const auto& entry : std::filesystem::directory_iterator(pathToDir)) {
  //   if (!entry.is_regular_file() || entry.path().extension() != ".root") {
  //     continue;
  //   }
  //   std::string pathToFile = entry.path();
  //   readerCfg.filePaths.push_back(pathToFile);
  // }

  // // Add the reader to the sequencer
  // sequencer.addReader(
  //     std::make_shared<RootSimClusterReader>(readerCfg, logLevel));

  // --------------------------------------------------------------
  // Ideal seeding setup

  IdealSeedingAlgorithm::Config idealSeedingCfg;
  idealSeedingCfg.inputSourceLinks = "Measurements";
  idealSeedingCfg.inputSimClusters = "SimClusters";
  idealSeedingCfg.outputSeeds = "IdealSeeds";
  idealSeedingCfg.minLayers = 10;
  idealSeedingCfg.maxLayers = 10;
  idealSeedingCfg.minSeedSize = 10;
  idealSeedingCfg.maxSeedSize = 10;

  sequencer.addAlgorithm(
      std::make_shared<IdealSeedingAlgorithm>(idealSeedingCfg, logLevel));

  // --------------------------------------------------------------
  // HT seeding setup

  HoughTransformSeeder::Config htSeederCfg;
  htSeederCfg.boundBoxHalfX = ag::GeometryOptions::instance()->tc1HalfPrimary;
  htSeederCfg.boundBoxHalfY = ag::GeometryOptions::instance()->tcHalfLong;
  htSeederCfg.boundBoxHalfZ = ag::GeometryOptions::instance()->tcHalfShort;

  htSeederCfg.nCellsX = 1000;
  htSeederCfg.nCellsY = 1000;

  htSeederCfg.minSeedSize = 5;
  htSeederCfg.maxSeedSize = 1000;

  htSeederCfg.nLSIterations = 20;

  ApollonSeedingAlgorithm::Config seedingAlgoCfg;
  seedingAlgoCfg.htSeeder = std::make_shared<HoughTransformSeeder>(htSeederCfg);
  seedingAlgoCfg.inputSourceLinks = "Measurements";
  seedingAlgoCfg.outputSeeds = "Seeds";
  seedingAlgoCfg.minLayers = 10;
  seedingAlgoCfg.maxLayers = 10;
  seedingAlgoCfg.minSeedSize = 10;
  seedingAlgoCfg.maxSeedSize = 100;
  seedingAlgoCfg.minScanEnergy = 0.2_GeV;
  seedingAlgoCfg.maxScanEnergy = 0.4_GeV;
  seedingAlgoCfg.energyScanStep = 0.01_GeV;
  seedingAlgoCfg.maxConnectionDistance = 1e16;

  sequencer.addAlgorithm(
      std::make_shared<ApollonSeedingAlgorithm>(seedingAlgoCfg, logLevel));

  // --------------------------------------------------------------
  // Track finding

  Acts::Experimental::DetectorNavigator::Config ckfNavigatorCfg;
  ckfNavigatorCfg.detector = detector.get();
  ckfNavigatorCfg.resolvePassive = false;
  ckfNavigatorCfg.resolveMaterial = true;
  ckfNavigatorCfg.resolveSensitive = true;
  Acts::Experimental::DetectorNavigator ckfNavigator(
      ckfNavigatorCfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));

  Acts::EigenStepper<> ckfStepper(field);
  auto ckfPropagator =
      Propagator(std::move(ckfStepper), std::move(ckfNavigator),
                 Acts::getDefaultLogger("Propagator", logLevel));

  Acts::CombinatorialKalmanFilter<Propagator, CKFTrackContainer> ckf(
      ckfPropagator,
      Acts::getDefaultLogger("CombinatorialKalmanFilter", logLevel));

  // Configuration for the measurement selector
  std::vector<
      std::pair<Acts::GeometryIdentifier, Acts::MeasurementSelectorCuts>>
      cuts;
  for (auto& vol : detector->volumes()) {
    for (auto& surf : vol->surfaces()) {
      if (!surf->geometryId().sensitive()) {
        continue;
      }
      cuts.push_back({surf->geometryId(),
                      {{},
                       {std::numeric_limits<double>::max()},
                       {1u},
                       {std::numeric_limits<double>::max()}}});
    }
  }
  Acts::MeasurementSelector::Config measurementSelectorCfg(cuts);

  Acts::MeasurementSelector measSel{measurementSelectorCfg};

  // CKF extensions
  Acts::GainMatrixUpdater ckfUpdater;

  Acts::CombinatorialKalmanFilterExtensions<CKFTrackContainer> ckfExtensions;
  ckfExtensions.calibrator.template connect<
      &simpleSourceLinkCalibrator<TrackStateContainerBackend>>();
  ckfExtensions.updater.template connect<
      &Acts::GainMatrixUpdater::operator()<TrackStateContainerBackend>>(
      &ckfUpdater);
  ckfExtensions.measurementSelector.template connect<
      &Acts::MeasurementSelector::select<TrackStateContainerBackend>>(&measSel);

  CKFTrackFindingAlgorithm::Config trackFindingCfg{
      .ckf = ckf,
  };
  trackFindingCfg.extensions = ckfExtensions;
  trackFindingCfg.inputSeeds = "Seeds";
  trackFindingCfg.outputTrackCandidates = "TrackCandidates";
  trackFindingCfg.outputTrackView = "CandidatesTrackView";
  trackFindingCfg.minCandidateSize = 10;
  trackFindingCfg.maxCandidateSize = 10;
  trackFindingCfg.maxSteps = 1e4;

  // auto trackFindingAlgorithm =
  //     std::make_shared<CKFTrackFindingAlgorithm>(trackFindingCfg, logLevel);
  // sequencer.addAlgorithm(trackFindingAlgorithm);

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
      Acts::AngleAxis3(ag::GeometryOptions::instance()->toWorldAngleX,
                       Acts::Vector3::UnitX())
          .toRotationMatrix();
  Acts::RotationMatrix3 refSurfToWorldRotationY =
      Acts::AngleAxis3(ag::GeometryOptions::instance()->toWorldAngleY,
                       Acts::Vector3::UnitY())
          .toRotationMatrix();
  Acts::RotationMatrix3 refSurfToWorldRotationZ =
      Acts::AngleAxis3(ag::GeometryOptions::instance()->toWorldAngleZ,
                       Acts::Vector3::UnitZ())
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
  KFTrackFittingAlgorithm::Config fitterCfg{
      .inputTrackCandidates = "TrackCandidates",
      .outputTracks = "Tracks",
      .fitter = fitter,
      .kfOptions = options};

  // sequencer.addAlgorithm(
  //     std::make_shared<KFTrackFittingAlgorithm>(fitterCfg, logLevel));

  // --------------------------------------------------------------
  // Event write out

  // Sim cluster writer
  auto clusterWriterCfg = RootSimClusterWriter::Config();
  clusterWriterCfg.inputClusters = "SimClusters";
  clusterWriterCfg.treeName = "clusters";
  clusterWriterCfg.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/test/clusters.root";

  sequencer.addWriter(
      std::make_shared<RootSimClusterWriter>(clusterWriterCfg, logLevel));

  // Seed writer
  auto seedWriterCfg = RootSimSeedWriter::Config();
  seedWriterCfg.inputSeeds = "Seeds";
  seedWriterCfg.inputTruthClusters = "SimClusters";
  seedWriterCfg.treeName = "seeds";
  seedWriterCfg.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/test/"
      "seeds.root";

  sequencer.addWriter(
      std::make_shared<RootSimSeedWriter>(seedWriterCfg, logLevel));

  auto idealSeedWriterCfg = RootSimSeedWriter::Config();
  idealSeedWriterCfg.inputSeeds = "IdealSeeds";
  idealSeedWriterCfg.inputTruthClusters = "SimClusters";
  idealSeedWriterCfg.treeName = "seeds";
  idealSeedWriterCfg.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/test/"
      "ideal-seeds.root";

  sequencer.addWriter(
      std::make_shared<RootSimSeedWriter>(idealSeedWriterCfg, logLevel));

  // Track candidate writer
  auto trackCandidateWriterCfg = RootSimTrackCandidateWriter::Config();
  trackCandidateWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  trackCandidateWriterCfg.inputTrackCandidates = "CandidatesTrackView";
  trackCandidateWriterCfg.treeName = "track-candidates";
  trackCandidateWriterCfg.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/test/"
      "track-candidates.root";

  // sequencer.addWriter(std::make_shared<RootSimTrackCandidateWriter>(
  //     trackCandidateWriterCfg, logLevel));

  // Fitted track writer
  auto trackWriterCfg = RootSimTrackWriter::Config();
  trackWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);
  trackWriterCfg.inputTracks = "Tracks";
  trackWriterCfg.inputTruthClusters = "SimClusters";
  trackWriterCfg.treeName = "fitted-tracks";
  trackWriterCfg.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/test/"
      "fitted-tracks.root";

  // sequencer.addWriter(
  //     std::make_shared<RootSimTrackWriter>(trackWriterCfg, logLevel));

  return sequencer.run();
}
