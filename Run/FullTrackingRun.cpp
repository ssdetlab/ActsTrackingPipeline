#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/GlobalChiSquareFitter.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <unordered_map>

#include <unistd.h>

#include "TrackingPipeline/Geometry/ApollonGeometry.hpp"
#include "TrackingPipeline/Geometry/ApollonGeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/ApollonRootSimDataReader.hpp"
#include "TrackingPipeline/Io/RootSeedWriter.hpp"
#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"
#include "TrackingPipeline/Io/RootTrackCandidateWriter.hpp"
#include "TrackingPipeline/Io/RootTrackWriter.hpp"
#include "TrackingPipeline/TrackFinding/CKFTrackFindingAlgorithm.hpp"
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
  Acts::Logging::Level logLevel = Acts::Logging::INFO;

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
  Sequencer sequencer(seqCfg);

  // Add the sim data reader
  ApollonIo::ApollonRootSimDataReader::Config readerCfg =
      ApollonIo::defaultSimConfig();
  readerCfg.outputSourceLinks = "Measurements";
  readerCfg.outputSimClusters = "SimClusters";
  readerCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);
  std::string pathToDir =
      "/home/romanurmanov/work/Apollon/geant4_sims/al_window_flange/out_data/"
      "temp";

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

  // --------------------------------------------------------------
  // Ideal seeding setup

  // sequencer.addAlgorithm(
  //     std::make_shared<PathSeedingAlgorithm>(seedingAlgoCfg, logLevel));

  // // --------------------------------------------------------------
  // // Track finding

  // Acts::Experimental::DetectorNavigator::Config ckfNavigatorCfg;
  // ckfNavigatorCfg.detector = detector.get();
  // ckfNavigatorCfg.resolvePassive = false;
  // ckfNavigatorCfg.resolveMaterial = true;
  // ckfNavigatorCfg.resolveSensitive = true;
  // Acts::Experimental::DetectorNavigator ckfNavigator(
  //     ckfNavigatorCfg, Acts::getDefaultLogger("DetectorNavigator",
  //     logLevel));

  // Acts::EigenStepper<> ckfStepper(field);
  // auto ckfPropagator =
  //     Propagator(std::move(ckfStepper), std::move(ckfNavigator),
  //                Acts::getDefaultLogger("Propagator", logLevel));

  // Acts::CombinatorialKalmanFilter<Propagator, CKFTrackContainer> ckf(
  //     ckfPropagator,
  //     Acts::getDefaultLogger("CombinatorialKalmanFilter", logLevel));

  // // Configuration for the measurement selector
  // std::vector<
  //     std::pair<Acts::GeometryIdentifier, Acts::MeasurementSelectorCuts>>
  //     cuts;
  // for (auto& vol : detector->volumes()) {
  //   for (auto& surf : vol->surfaces()) {
  //     if (vol->name() == "layer0") {
  //       cuts.push_back({surf->geometryId(),
  //                       {{}, {std::numeric_limits<double>::max()},
  //                       {1000u}}});
  //     } else {
  //       cuts.push_back({surf->geometryId(),
  //                       {{},
  //                        {std::numeric_limits<double>::max()},
  //                        {1u},
  //                        {std::numeric_limits<double>::max()}}});
  //     }
  //   }
  // }
  // Acts::MeasurementSelector::Config measurementSelectorCfg(cuts);

  // Acts::MeasurementSelector measSel{measurementSelectorCfg};

  // // CKF extensions
  // Acts::GainMatrixUpdater ckfUpdater;

  // Acts::CombinatorialKalmanFilterExtensions<CKFTrackContainer> ckfExtensions;
  // ckfExtensions.calibrator.template connect<
  //     &simpleSourceLinkCalibrator<TrackStateContainerBackend>>();
  // ckfExtensions.updater.template connect<
  //     &Acts::GainMatrixUpdater::operator()<TrackStateContainerBackend>>(
  //     &ckfUpdater);
  // ckfExtensions.measurementSelector.template connect<
  //     &Acts::MeasurementSelector::select<TrackStateContainerBackend>>(&measSel);

  // CKFTrackFindingAlgorithm::Config trackFindingCfg{
  //     .ckf = ckf,
  // };
  // trackFindingCfg.extensions = ckfExtensions;
  // trackFindingCfg.inputSeeds = "PathSeeds";
  // trackFindingCfg.outputTrackCandidates = "TrackCandidates";
  // trackFindingCfg.outputTrackView = "CandidatesTrackView";
  // trackFindingCfg.minCandidateSize = 5;
  // trackFindingCfg.maxCandidateSize = 5;
  // trackFindingCfg.maxSteps = 1e4;

  // auto trackFindingAlgorithm =
  //     std::make_shared<CKFTrackFindingAlgorithm>(trackFindingCfg, logLevel);
  // sequencer.addAlgorithm(trackFindingAlgorithm);

  // // --------------------------------------------------------------
  // // Track fitting

  // Acts::GainMatrixUpdater kfUpdater;
  // Acts::GainMatrixSmoother kfSmoother;

  // // Initialize track fitter options
  // Acts::KalmanFitterExtensions<RecoTrajectory> extensions;
  // // Add calibrator
  // extensions.calibrator.connect<&simpleSourceLinkCalibrator<RecoTrajectory>>();
  // // Add the updater
  // extensions.updater
  //     .connect<&Acts::GainMatrixUpdater::operator()<RecoTrajectory>>(
  //         &kfUpdater);
  // // Add the smoother
  // extensions.smoother
  //     .connect<&Acts::GainMatrixSmoother::operator()<RecoTrajectory>>(
  //         &kfSmoother);
  // // Add the surface accessor
  // extensions.surfaceAccessor
  //     .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
  //         &surfaceAccessor);

  // auto propOptions = PropagatorOptions(gctx, mctx);

  // propOptions.maxSteps = 1e5;

  // auto options =
  //     Acts::KalmanFitterOptions(gctx, mctx, cctx, extensions, propOptions);

  // // Reference surface for sampling the track
  // double halfX = std::numeric_limits<double>::max();
  // double halfY = std::numeric_limits<double>::max();

  // // double refZ = 13060.6_mm + 914_mm / 2;
  // double refZ = gOpt.beWindowTranslation[2];
  // // double refZ = gOpt.staveZ.at(8) - 0.1_mm;
  // Acts::Transform3 transform(Acts::Translation3(Acts::Vector3(0, refZ, 0)) *
  //                            gOpt.actsToWorldRotation.inverse());

  // auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
  //     transform, std::make_shared<Acts::RectangleBounds>(halfX, halfY));

  // Acts::GeometryIdentifier geoId;
  // geoId.setExtra(1);
  // refSurface->assignGeometryId(std::move(geoId));

  // options.referenceSurface = refSurface.get();

  // Acts::Experimental::DetectorNavigator::Config cfg;
  // cfg.detector = detector.get();
  // cfg.resolvePassive = false;
  // cfg.resolveMaterial = true;
  // cfg.resolveSensitive = true;
  // Acts::Experimental::DetectorNavigator kfNavigator(
  //     cfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));

  // Acts::EigenStepper<> kfStepper(std::move(field));
  // auto kfPropagator =
  //     Propagator(std::move(kfStepper), std::move(kfNavigator),
  //                Acts::getDefaultLogger("Propagator", logLevel));

  // const auto fitter = KF(
  //     kfPropagator, Acts::getDefaultLogger("DetectorKalmanFilter",
  //     logLevel));

  // // Add the track fitting algorithm to the sequencer
  // KFTrackFittingAlgorithm::Config fitterCfg{
  //     .inputTrackCandidates = "TrackCandidates",
  //     .outputTracks = "Tracks",
  //     .fitter = fitter,
  //     .kfOptions = options};

  // sequencer.addAlgorithm(
  //     std::make_shared<KFTrackFittingAlgorithm>(fitterCfg, logLevel));

  // --------------------------------------------------------------
  // Event write out

  // --------------------------------------------------------------
  // Event write out

  // Sim cluster writer
  auto clusterWriterCfg = RootSimClusterWriter::Config();

  clusterWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  clusterWriterCfg.inputClusters = "Clusters";
  clusterWriterCfg.treeName = "clusters";
  clusterWriterCfg.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/test/clusters.root";

  // sequencer.addWriter(
  //     std::make_shared<RootSimClusterWriter>(clusterWriterCfg, logLevel));

  // Seed writer
  auto seedWriterCfg = RootSeedWriter::Config();

  seedWriterCfg.inputSeeds = "PathSeeds";
  seedWriterCfg.treeName = "seeds";
  seedWriterCfg.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/test/"
      "seeds-data.root";

  seedWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  // sequencer.addWriter(
  //     std::make_shared<RootSeedWriter>(seedWriterCfg, logLevel));

  // Track candidate writer
  auto trackCandidateWriterCfg = RootTrackCandidateWriter::Config();
  trackCandidateWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  trackCandidateWriterCfg.inputTrackCandidates = "CandidatesTrackView";
  trackCandidateWriterCfg.treeName = "track-candidates";
  trackCandidateWriterCfg.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/test/"
      "track-candidates-data.root";

  // sequencer.addWriter(std::make_shared<RootTrackCandidateWriter>(
  //     trackCandidateWriterCfg, logLevel));

  // Fitted track writer
  auto trackWriterCfg = RootTrackWriter::Config();
  trackWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  trackWriterCfg.inputTracks = "Tracks";
  trackWriterCfg.treeName = "fitted-tracks";
  trackWriterCfg.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/test/"
      "fitted-tracks-data.root";

  // sequencer.addWriter(
  //     std::make_shared<RootTrackWriter>(trackWriterCfg, logLevel));

  return sequencer.run();
}
