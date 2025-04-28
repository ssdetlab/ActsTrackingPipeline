#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Seeding/PathSeeder.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <filesystem>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/E320RootDataReader.hpp"
#include "TrackingPipeline/Io/JsonTrackLookupReader.hpp"
#include "TrackingPipeline/Io/RootMeasurementWriter.hpp"
#include "TrackingPipeline/Io/RootSimSeedWriter.hpp"
#include "TrackingPipeline/Io/RootSimTrackCandidateWriter.hpp"
#include "TrackingPipeline/Io/RootSimTrackWriter.hpp"
#include "TrackingPipeline/Io/RootTrackParamsReader.hpp"
#include "TrackingPipeline/MagneticField/CompositeMagField.hpp"
#include "TrackingPipeline/MagneticField/ConstantBoundedField.hpp"
#include "TrackingPipeline/MagneticField/DipoleMagField.hpp"
#include "TrackingPipeline/MagneticField/QuadrupoleMagField.hpp"
#include "TrackingPipeline/TrackFinding/CKFTrackFindingAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/ForwardOrderedIntersectionFinder.hpp"
#include "TrackingPipeline/TrackFinding/LayerPathWidthProvider.hpp"
#include "TrackingPipeline/TrackFinding/PathSeedingAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/E320SourceLinkGridConstructor.hpp"
#include "TrackingPipeline/TrackFinding/TrackLookupProvider.hpp"
#include "TrackingPipeline/TrackFitting/TrackFittingAlgorithm.hpp"

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
using KFTrajectory = TrackFittingAlgorithm::Trajectory;
using KFTrackContainer = TrackFittingAlgorithm::TrackContainer;
using KF = Acts::KalmanFitter<Propagator, KFTrajectory>;

using namespace Acts::UnitLiterals;

int main() {
  // Set the log level
  Acts::Logging::Level logLevel = Acts::Logging::VERBOSE;

  // Dummy context and options
  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext cctx;
  E320Geometry::GeometryOptions gOpt;

  // --------------------------------------------------------------
  // Detector setup

  // Set the path to the gdml file
  // and the names of the volumes to be converted
  std::string gdmlPath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_gdmls/"
      "ett_geometry_f566a577.gdml";
  std::vector<std::string> names{"OPPPSensitive", "DetChamberWindow"};

  // Veto PDC window material mapping
  // to preserve homogeneous material
  // from Geant4
  Acts::GeometryIdentifier pdcWindowId;
  pdcWindowId.setApproach(1);
  std::vector<Acts::GeometryIdentifier> materialVeto{pdcWindowId};

  std::string materialPath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_material/"
      "Uniform_DirectZ_TrackerOnly_256x128_1M/material.json";

  // Build the detector
  auto trackerBP = E320Geometry::makeBlueprintE320(gdmlPath, names, gOpt);
  auto detector = E320Geometry::buildE320Detector(
      std::move(trackerBP), gctx, gOpt, materialPath, materialVeto);

  // --------------------------------------------------------------
  // The magnetic field setup

  // Extent in already rotated frame
  Acts::Extent quad1Extent;
  quad1Extent.set(Acts::BinningValue::binX,
                  gOpt.quad1Translation[0] - gOpt.quad1Bounds[0],
                  gOpt.quad1Translation[0] + gOpt.quad1Bounds[0]);
  quad1Extent.set(Acts::BinningValue::binZ,
                  gOpt.quad1Translation[1] - gOpt.quad1Bounds[1],
                  gOpt.quad1Translation[1] + gOpt.quad1Bounds[1]);
  quad1Extent.set(Acts::BinningValue::binY,
                  gOpt.quad1Translation[2] - gOpt.quad1Bounds[2],
                  gOpt.quad1Translation[2] + gOpt.quad1Bounds[2]);

  Acts::Extent quad2Extent;
  quad2Extent.set(Acts::BinningValue::binX,
                  gOpt.quad2Translation[0] - gOpt.quad2Bounds[0],
                  gOpt.quad2Translation[0] + gOpt.quad2Bounds[0]);
  quad2Extent.set(Acts::BinningValue::binZ,
                  gOpt.quad2Translation[1] - gOpt.quad2Bounds[1],
                  gOpt.quad2Translation[1] + gOpt.quad2Bounds[1]);
  quad2Extent.set(Acts::BinningValue::binY,
                  gOpt.quad2Translation[2] - gOpt.quad2Bounds[2],
                  gOpt.quad2Translation[2] + gOpt.quad2Bounds[2]);

  Acts::Extent quad3Extent;
  quad3Extent.set(Acts::BinningValue::binX,
                  gOpt.quad3Translation[0] - gOpt.quad3Bounds[0],
                  gOpt.quad3Translation[0] + gOpt.quad3Bounds[0]);
  quad3Extent.set(Acts::BinningValue::binZ,
                  gOpt.quad3Translation[1] - gOpt.quad3Bounds[1],
                  gOpt.quad3Translation[1] + gOpt.quad3Bounds[1]);
  quad3Extent.set(Acts::BinningValue::binY,
                  gOpt.quad3Translation[2] - gOpt.quad3Bounds[2],
                  gOpt.quad3Translation[2] + gOpt.quad3Bounds[2]);

  Acts::Extent dipoleExtent;
  dipoleExtent.set(Acts::BinningValue::binX,
                   gOpt.dipoleTranslation.x() - gOpt.dipoleBounds[0],
                   gOpt.dipoleTranslation.x() + gOpt.dipoleBounds[0]);
  dipoleExtent.set(Acts::BinningValue::binZ,
                   gOpt.dipoleTranslation.y() - gOpt.dipoleBounds[1],
                   gOpt.dipoleTranslation.y() + gOpt.dipoleBounds[1]);
  dipoleExtent.set(Acts::BinningValue::binY,
                   gOpt.dipoleTranslation.z() - gOpt.dipoleBounds[2],
                   gOpt.dipoleTranslation.z() + gOpt.dipoleBounds[2]);

  Acts::Extent xCorrectorExtent;
  xCorrectorExtent.set(
      Acts::BinningValue::binX,
      gOpt.xCorrectorTranslation.x() - gOpt.xCorrectorBounds[0],
      gOpt.xCorrectorTranslation.x() + gOpt.xCorrectorBounds[0]);
  xCorrectorExtent.set(
      Acts::BinningValue::binZ,
      gOpt.xCorrectorTranslation.y() - gOpt.xCorrectorBounds[1],
      gOpt.xCorrectorTranslation.y() + gOpt.xCorrectorBounds[1]);
  xCorrectorExtent.set(
      Acts::BinningValue::binY,
      gOpt.xCorrectorTranslation.z() - gOpt.xCorrectorBounds[2],
      gOpt.xCorrectorTranslation.z() + gOpt.xCorrectorBounds[2]);

  QuadrupoleMagField quad1Field(
      gOpt.quadrupolesParams[0],
      gOpt.actsToWorldRotation.inverse() * gOpt.quad1Translation,
      gOpt.actsToWorldRotation);
  QuadrupoleMagField quad2Field(
      gOpt.quadrupolesParams[1],
      gOpt.actsToWorldRotation.inverse() * gOpt.quad2Translation,
      gOpt.actsToWorldRotation);
  QuadrupoleMagField quad3Field(
      gOpt.quadrupolesParams[2],
      gOpt.actsToWorldRotation.inverse() * gOpt.quad3Translation,
      gOpt.actsToWorldRotation);

  double dipoleB = 0.2192_T;
  DipoleMagField dipoleField(
      gOpt.dipoleParams, dipoleB, gOpt.actsToWorldRotation,
      gOpt.actsToWorldRotation.inverse() * gOpt.dipoleTranslation);

  // TODO: Add the real field value
  Acts::Vector3 xCorrectorB(0, 0, 0.32_T);
  ConstantBoundedField xCorrectorField(xCorrectorB, xCorrectorExtent);

  CompositeMagField::FieldComponents fieldComponents = {
      {quad1Extent, &quad1Field},
      {quad2Extent, &quad2Field},
      {quad3Extent, &quad3Field},
      {dipoleExtent, &dipoleField},
      {xCorrectorExtent, &xCorrectorField}};

  auto field = std::make_shared<CompositeMagField>(fieldComponents);

  // --------------------------------------------------------------
  // Event reading
  SimpleSourceLink::SurfaceAccessor surfaceAccessor{detector.get()};

  // Setup the sequencer
  Sequencer::Config seqCfg;
  seqCfg.events = 1e3;
  seqCfg.numThreads = 1;
  seqCfg.trackFpes = false;
  Sequencer sequencer(seqCfg);

  // Add the sim data reader
  E320Io::E320RootDataReader::Config readerCfg;
  //   readerCfg.clusterFilter = hourglassFilter;
  readerCfg.treeName = "MyTree";
  readerCfg.outputSourceLinks = "Measurements";
  std::string pathToDir =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_dataInRootFormat/"
      "E320Shift_Nov_2024";

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
      std::make_shared<E320Io::E320RootDataReader>(readerCfg, logLevel));

  //  // --------------------------------------------------------------
  //  // The path seeding setup
  //  auto pathSeederCfg = Acts::PathSeeder::Config();
  //
  //  // Combine layers into surfaces to take care of the gaps
  //  std::vector<std::shared_ptr<Acts::Surface>> layerPtrs;
  //  for (int i = 0; i < gOpt.staveZ.size(); i++) {
  //    // Bounds
  //    double halfX = gOpt.chipSizeX / 2;
  //    double halfY = ((gOpt.chipY.at(8) + gOpt.chipSizeY / 2) -
  //                    (gOpt.chipY.at(0) - gOpt.chipSizeY / 2)) /
  //                   2;
  //
  //    // Center
  //    double centerY = ((-gOpt.chipY.at(8) - gOpt.chipSizeY / 2) +
  //                      (-gOpt.chipY.at(0) + gOpt.chipSizeY / 2)) /
  //                     2;
  //
  //    Acts::Transform3 transform(Acts::Translation3(Acts::Vector3(
  //                                   gOpt.chipX, gOpt.staveZ.at(i), centerY))
  //                                   *
  //                               gOpt.actsToWorldRotation.inverse());
  //
  //    auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
  //        transform, std::make_shared<Acts::RectangleBounds>(halfX, halfY));
  //
  //    Acts::GeometryIdentifier geoId;
  //    geoId.setSensitive(i + 1);
  //    surface->assignGeometryId(std::move(geoId));
  //    layerPtrs.push_back(surface);
  //
  //    if (i == 0) {
  //      pathSeederCfg.refLayerIds.push_back(geoId);
  //    }
  //  }
  //
  //  // Intersection finder
  //  std::vector<const Acts::Surface*> layers;
  //  for (auto layer : layerPtrs) {
  //    layers.push_back(layer.get());
  //  }
  //
  //  ForwardOrderedIntersectionFinder::Config intersectionFinderCfg;
  //  intersectionFinderCfg.layers = std::move(layers);
  //  intersectionFinderCfg.tol = (gOpt.chipY.at(1) - gOpt.chipSizeY / 2) -
  //                              (gOpt.chipY.at(0) + gOpt.chipSizeY / 2) +
  //                              1_mm;
  //  ForwardOrderedIntersectionFinder
  //  intersectionFinder(intersectionFinderCfg);
  //
  //  pathSeederCfg.intersectionFinder
  //      .connect<&ForwardOrderedIntersectionFinder::operator()>(
  //          &intersectionFinder);
  //
  //  // Path width provider
  //  std::map<std::int32_t, std::pair<double, double>> pathWidths = {
  //      {0, {100_um, 100_um}},
  //      {1, {300_um, 300_um}},
  //      {2, {305_um, 350_um}},
  //      {3, {310_um, 400_um}},
  //  };
  //
  //  LayerPathWidthProvider pathWidthProvider(pathWidths);
  //
  //  pathSeederCfg.pathWidthProvider.connect<&LayerPathWidthProvider::operator()>(
  //      &pathWidthProvider);
  //
  //  // Grid to bin the source links
  //  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
  //  layerMap; for (auto layer : layerPtrs) {
  //    layerMap.insert({layer->geometryId(), layer.get()});
  //  }
  //
  //  SourceLinkGridConstructor::Config gridConstructorCfg{
  //      .bins = std::make_pair(1, 1000), .layers = layerMap};
  //  gridConstructorCfg.surfaceAccessor
  //      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
  //          &surfaceAccessor);
  //
  //  auto gridConstructor =
  //      std::make_shared<SourceLinkGridConstructor>(gridConstructorCfg);
  //
  //  // Estimator of the IP and first hit
  //  // parameters of the track
  //  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
  //  refLayers; const auto& refVolume = detector->findDetectorVolume("layer0");
  //  for (const auto* surf : refVolume->surfaces()) {
  //    refLayers.try_emplace(surf->geometryId(), surf);
  //  }
  //
  //  JsonTrackLookupReader::Config lookupReaderCfg;
  //  lookupReaderCfg.refLayers = refLayers;
  //  lookupReaderCfg.bins = {1000, 1};
  //
  //  TrackLookupProvider::Config lookupProviderCfg;
  //  lookupProviderCfg.lookupPath =
  //      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_lookups/"
  //      "RangedUniform_05_45_Stationary_000_1000x1_200k/lookup.json";
  //  lookupProviderCfg.trackLookupReader =
  //      std::make_shared<JsonTrackLookupReader>(lookupReaderCfg);
  //  TrackLookupProvider lookupProvider(lookupProviderCfg);
  //
  //  pathSeederCfg.trackEstimator.connect<&TrackLookupProvider::lookup>(
  //      &lookupProvider);
  //
  //  // Create the path seeder algorithm
  //  auto seedingAlgoCfg = PathSeedingAlgorithm::Config();
  //  seedingAlgoCfg.seeder = std::make_shared<Acts::PathSeeder>(pathSeederCfg);
  //  seedingAlgoCfg.sourceLinkGridConstructor = gridConstructor;
  //  seedingAlgoCfg.inputSourceLinks = "Measurements";
  //  seedingAlgoCfg.outputSeeds = "PathSeeds";
  //  seedingAlgoCfg.minSeedSize = 4;
  //  seedingAlgoCfg.maxSeedSize = 100;
  //  sequencer.addAlgorithm(
  //      std::make_shared<PathSeedingAlgorithm>(seedingAlgoCfg, logLevel));
  //
  //  // --------------------------------------------------------------
  //  // Track finding
  //  Acts::Experimental::DetectorNavigator::Config ckfNavigatorCfg;
  //  ckfNavigatorCfg.detector = detector.get();
  //  ckfNavigatorCfg.resolvePassive = false;
  //  ckfNavigatorCfg.resolveMaterial = true;
  //  ckfNavigatorCfg.resolveSensitive = true;
  //  Acts::Experimental::DetectorNavigator ckfNavigator(
  //      ckfNavigatorCfg, Acts::getDefaultLogger("DetectorNavigator",
  //      logLevel));
  //
  //  Acts::EigenStepper<> ckfStepper(field);
  //  auto ckfPropagator =
  //      Propagator(std::move(ckfStepper), std::move(ckfNavigator),
  //                 Acts::getDefaultLogger("Propagator", logLevel));
  //
  //  Acts::CombinatorialKalmanFilter<Propagator, CKFTrackContainer> ckf(
  //      ckfPropagator,
  //      Acts::getDefaultLogger("CombinatorialKalmanFilter", logLevel));
  //
  //  // Configuration for the measurement selector
  //  std::vector<
  //      std::pair<Acts::GeometryIdentifier, Acts::MeasurementSelectorCuts>>
  //      cuts;
  //  for (auto& vol : detector->volumes()) {
  //    for (auto& surf : vol->surfaces()) {
  //      if (vol->name() == "layer0") {
  //        cuts.push_back({surf->geometryId(),
  //                        {{}, {std::numeric_limits<double>::max()},
  //                        {1000u}}});
  //      } else {
  //        cuts.push_back({surf->geometryId(), {{}, {15}, {1u}}});
  //      }
  //    }
  //  }
  //  Acts::MeasurementSelector::Config measurementSelectorCfg(cuts);
  //
  //  Acts::MeasurementSelector measSel{measurementSelectorCfg};
  //
  //  // CKF extensions
  //  Acts::GainMatrixUpdater ckfUpdater;
  //
  //  Acts::CombinatorialKalmanFilterExtensions<CKFTrackContainer>
  //  ckfExtensions; ckfExtensions.calibrator.template connect<
  //      &simpleSourceLinkCalibrator<TrackStateContainerBackend>>();
  //  ckfExtensions.updater.template connect<
  //      &Acts::GainMatrixUpdater::operator()<TrackStateContainerBackend>>(
  //      &ckfUpdater);
  //  ckfExtensions.measurementSelector.template connect<
  //      &Acts::MeasurementSelector::select<TrackStateContainerBackend>>(&measSel);
  //
  //  CKFTrackFindingAlgorithm::Config trackFindingCfg{
  //      .ckf = ckf,
  //  };
  //  trackFindingCfg.extensions = ckfExtensions;
  //  trackFindingCfg.inputSeeds = "PathSeeds";
  //  trackFindingCfg.outputTrackCandidates = "TrackCandidates";
  //  trackFindingCfg.outputTrackView = "CandidatesTrackView";
  //  trackFindingCfg.minCandidateSize = 4;
  //  trackFindingCfg.maxCandidateSize = 4;
  //  trackFindingCfg.maxSteps = 1e4;
  //
  //  auto trackFindingAlgorithm =
  //      std::make_shared<CKFTrackFindingAlgorithm>(trackFindingCfg, logLevel);
  //  sequencer.addAlgorithm(trackFindingAlgorithm);
  //
  //  // --------------------------------------------------------------
  //  // Track fitting
  //
  //  Acts::GainMatrixUpdater kfUpdater;
  //  Acts::GainMatrixSmoother kfSmoother;
  //
  //  // Initialize track fitter options
  //  Acts::KalmanFitterExtensions<KFTrajectory> extensions;
  //  // Add calibrator
  //  extensions.calibrator.connect<&simpleSourceLinkCalibrator<KFTrajectory>>();
  //  // Add the updater
  //  extensions.updater
  //      .connect<&Acts::GainMatrixUpdater::operator()<KFTrajectory>>(&kfUpdater);
  //  // Add the smoother
  //  extensions.smoother
  //      .connect<&Acts::GainMatrixSmoother::operator()<KFTrajectory>>(
  //          &kfSmoother);
  //  // Add the surface accessor
  //  extensions.surfaceAccessor
  //      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
  //          &surfaceAccessor);
  //
  //  auto propOptions = PropagatorOptions(gctx, mctx);
  //
  //  propOptions.maxSteps = 300;
  //
  //  auto options =
  //      Acts::KalmanFitterOptions(gctx, mctx, cctx, extensions, propOptions);
  //
  //  // Reference surface for sampling the track at the IP
  //  double halfX = std::numeric_limits<double>::max();
  //  double halfY = std::numeric_limits<double>::max();
  //
  //  // double refZ = gOpt.dipoleTranslation.z() + gOpt.dipoleBounds.at(2);
  //  double refZ = 0;
  //  Acts::Transform3 transform(Acts::Translation3(Acts::Vector3(0, refZ, 0)) *
  //                             gOpt.actsToWorldRotation.inverse());
  //
  //  auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
  //      transform, std::make_shared<Acts::RectangleBounds>(halfX, halfY));
  //
  //  Acts::GeometryIdentifier geoId;
  //  geoId.setExtra(1);
  //  refSurface->assignGeometryId(std::move(geoId));
  //
  //  options.referenceSurface = refSurface.get();
  //
  //  Acts::Experimental::DetectorNavigator::Config cfg;
  //  cfg.detector = detector.get();
  //  cfg.resolvePassive = false;
  //  cfg.resolveMaterial = true;
  //  cfg.resolveSensitive = true;
  //  Acts::Experimental::DetectorNavigator kfNavigator(
  //      cfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));
  //
  //  Acts::EigenStepper<> kfStepper(std::move(field));
  //  auto kfPropagator =
  //      Propagator(std::move(kfStepper), std::move(kfNavigator),
  //                 Acts::getDefaultLogger("Propagator", logLevel));
  //
  //  const auto fitter = KF(
  //      kfPropagator, Acts::getDefaultLogger("DetectorKalmanFilter",
  //      logLevel));
  //
  //  // Add the track fitting algorithm to the sequencer
  //  TrackFittingAlgorithm::Config fitterCfg{
  //      .inputTrackCandidates = "TrackCandidates",
  //      .outputTracks = "Tracks",
  //      .fitter = fitter,
  //      .kfOptions = options};
  //
  //  sequencer.addAlgorithm(
  //      std::make_shared<TrackFittingAlgorithm>(fitterCfg, logLevel));
  //
  // --------------------------------------------------------------
  // Event write out

  // Sim cluster writer
  auto clusterWriterCfg = RootMeasurementWriter::Config();

  clusterWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  clusterWriterCfg.inputMeasurements = "Measurements";
  clusterWriterCfg.treeName = "clusters";
  clusterWriterCfg.filePath = "clusters.root";

  sequencer.addWriter(
      std::make_shared<RootMeasurementWriter>(clusterWriterCfg, logLevel));
  //
  //  // Seed writer
  //  auto seedWriterCfg = RootSimSeedWriter::Config();
  //
  //  seedWriterCfg.inputSeeds = "PathSeeds";
  //  seedWriterCfg.inputTruthClusters = "Clusters";
  //  seedWriterCfg.treeName = "seeds";
  //  seedWriterCfg.filePath = "seeds-sig-bkg.root";
  //  seedWriterCfg.targetTrueTrackSize = 4;
  //
  //  sequencer.addWriter(
  //      std::make_shared<RootSimSeedWriter>(seedWriterCfg, logLevel));
  //
  //  // Track candidate writer
  //  auto trackCandidateWriterCfg = RootSimTrackCandidateWriter::Config();
  //  trackCandidateWriterCfg.surfaceAccessor
  //      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
  //          &surfaceAccessor);
  //
  //  trackCandidateWriterCfg.inputTrackCandidates = "CandidatesTrackView";
  //  trackCandidateWriterCfg.inputTruthClusters = "Clusters";
  //  trackCandidateWriterCfg.treeName = "track-candidates";
  //  trackCandidateWriterCfg.filePath = "track-candidates-sig-bkg-test.root";
  //  trackCandidateWriterCfg.targetTrueTrackSize = 4;
  //
  //  sequencer.addWriter(std::make_shared<RootSimTrackCandidateWriter>(
  //      trackCandidateWriterCfg, logLevel));
  //
  //  // Fitted track writer
  //  auto trackWriterCfg = RootSimTrackWriter::Config();
  //  trackWriterCfg.surfaceAccessor
  //      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
  //          &surfaceAccessor);
  //
  //  trackWriterCfg.inputTracks = "Tracks";
  //  trackWriterCfg.inputTruthClusters = "Clusters";
  //  trackWriterCfg.treeName = "fitted-tracks";
  //  trackWriterCfg.filePath = "fitted-tracks-sig-bkg-test.root";
  //  trackWriterCfg.targetTrueTrackSize = 4;
  //
  //  sequencer.addWriter(
  //      std::make_shared<RootSimTrackWriter>(trackWriterCfg, logLevel));
  //
  return sequencer.run();
}
