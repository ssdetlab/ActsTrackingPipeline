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

#include <iostream>
#include <limits>
#include <memory>
#include <unordered_map>

#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"
#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/AlignmentResultWriter.hpp"
#include "TrackingPipeline/Io/DummyReader.hpp"
#include "TrackingPipeline/Io/JsonTrackLookupReader.hpp"
#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"
#include "TrackingPipeline/Io/RootSimSeedWriter.hpp"
#include "TrackingPipeline/Io/RootSimTrackCandidateWriter.hpp"
#include "TrackingPipeline/Io/RootSimTrackWriter.hpp"
#include "TrackingPipeline/Io/RootTrackParamsReader.hpp"
#include "TrackingPipeline/MagneticField/CompositeMagField.hpp"
#include "TrackingPipeline/MagneticField/ConstantBoundedField.hpp"
#include "TrackingPipeline/MagneticField/DipoleMagField.hpp"
#include "TrackingPipeline/MagneticField/QuadrupoleMagField.hpp"
#include "TrackingPipeline/Simulation/GaussianVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"
#include "TrackingPipeline/Simulation/MisalignedDigitizer.hpp"
#include "TrackingPipeline/Simulation/RangedUniformMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/StationaryVertexGenerator.hpp"
#include "TrackingPipeline/TrackFinding/CKFTrackFindingAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/E320SourceLinkGridConstructor.hpp"
#include "TrackingPipeline/TrackFinding/ForwardOrderedIntersectionFinder.hpp"
#include "TrackingPipeline/TrackFinding/LayerPathWidthProvider.hpp"
#include "TrackingPipeline/TrackFinding/PathSeedingAlgorithm.hpp"
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
  Acts::Logging::Level logLevel = Acts::Logging::INFO;

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

  auto aStore =
      std::make_shared<std::map<Acts::GeometryIdentifier, Acts::Transform3>>();
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive()) {
        aStore->emplace(s->geometryId(), s->transform(gctx));
      }
    }
  }
  AlignmentContext alignCtx(aStore);
  gctx = Acts::GeometryContext{alignCtx};

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
  Acts::Vector3 xCorrectorB(0, 0, 0);
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
  /*seqCfg.events = 1e3;*/
  seqCfg.numThreads = 1;
  seqCfg.trackFpes = false;
  Sequencer sequencer(seqCfg);

  // --------------------------------------------------------------
  // Add dummy reader
  DummyReader::Config dummyReaderCfg;
  dummyReaderCfg.outputSourceLinks = "SimMeasurements";
  dummyReaderCfg.outputSimClusters = "SimClusters";
  dummyReaderCfg.nEvents = 1;

  sequencer.addReader(std::make_shared<DummyReader>(dummyReaderCfg));

  // --------------------------------------------------------------
  // Setup the measurements creator

  Acts::Experimental::DetectorNavigator::Config navCfg;
  navCfg.detector = detector.get();
  navCfg.resolvePassive = false;
  navCfg.resolveMaterial = true;
  navCfg.resolveSensitive = true;

  Acts::Experimental::DetectorNavigator navigator(
      navCfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));
  Acts::EigenStepper<> stepper(field);

  auto propagator = Propagator(std::move(stepper), std::move(navigator));

  auto momGen = std::make_shared<RangedUniformMomentumGenerator>();
  momGen->Pranges = {
      {1.9_GeV, 2.1_GeV}, {2.1_GeV, 2.3_GeV}, {2.3_GeV, 2.5_GeV}};

  // Digitizer
  auto digitizer = std::make_shared<MisalignedDigitizer>();
  digitizer->resolution = {5_um, 5_um};
  digitizer->shifts = {{8, {0, 0}},
                       {6, {-12_um, 15_um}},
                       {4, {18_um, -15_um}},
                       {2, {-24_um, 20_um}},
                       {0, {0, 0}}};

  /*auto vertexGen = std::make_shared<GaussianVertexGenerator>(*/
  /*    Acts::Vector3(0.61872_mm, 16674_mm, -93.775_mm),*/
  /*    Acts::SquareMatrix3::Identity());*/

  auto vertexGen = std::make_shared<StationaryVertexGenerator>();
  vertexGen->vertex = Acts::Vector3(0, gOpt.beWindowTranslation[2], 0);

  MeasurementsCreator::Config mcCfg;
  mcCfg.vertexGenerator = vertexGen;
  mcCfg.momentumGenerator = momGen;
  mcCfg.hitDigitizer = digitizer;
  mcCfg.maxSteps = 300;
  mcCfg.isSignal = true;

  auto measurementsCreator =
      std::make_shared<MeasurementsCreator>(propagator, mcCfg);

  MeasurementsEmbeddingAlgorithm::Config mcaCfg;
  mcaCfg.inputSourceLinks = "SimMeasurements";
  mcaCfg.inputSimClusters = "SimClusters";
  mcaCfg.outputSourceLinks = "Measurements";
  mcaCfg.outputSimClusters = "Clusters";
  mcaCfg.measurementGenerator = measurementsCreator;
  mcaCfg.randomNumberSvc =
      std::make_shared<RandomNumbers>(RandomNumbers::Config());
  mcaCfg.nMeasurements = 5e2;

  sequencer.addAlgorithm(
      std::make_shared<MeasurementsEmbeddingAlgorithm>(mcaCfg, logLevel));

  // --------------------------------------------------------------
  // Reference surface for sampling the track at the IP
  double halfX = std::numeric_limits<double>::max();
  double halfY = std::numeric_limits<double>::max();

  double refZ = gOpt.beWindowTranslation[2];
  /*Acts::Transform3 transform(Acts::Translation3(Acts::Vector3(0, refZ, 0)) **/
  /*                           gOpt.actsToWorldRotation.inverse());*/
  Acts::Transform3 transform(
      Acts::Translation3(Acts::Vector3(0.61872_mm, 16674_mm, -93.775_mm)) *
      gOpt.actsToWorldRotation.inverse());

  auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      transform, std::make_shared<Acts::RectangleBounds>(halfX, halfY));
  Acts::GeometryIdentifier geoId;
  geoId.setExtra(1);
  refSurface->assignGeometryId(std::move(geoId));

  // --------------------------------------------------------------
  // The path seeding setup
  auto pathSeederCfg = Acts::PathSeeder::Config();

  // Intersection finder
  std::vector<const Acts::Surface*> layers;
  // Grid to bin the source links
  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> layerMap;
  for (auto& vol : detector->volumes()) {
    for (auto& surf : vol->surfaces()) {
      if (!surf->geometryId().sensitive()) {
        continue;
      }
      layers.push_back(surf);
      layerMap.insert({surf->geometryId(), surf});
    }
  }

  ForwardOrderedIntersectionFinder::Config intersectionFinderCfg;
  intersectionFinderCfg.layers = std::move(layers);
  ForwardOrderedIntersectionFinder intersectionFinder(intersectionFinderCfg);

  pathSeederCfg.intersectionFinder
      .connect<&ForwardOrderedIntersectionFinder::operator()>(
          &intersectionFinder);

  // Path width provider
  std::map<std::int32_t, std::pair<double, double>> pathWidths = {
      {8, {100_m, 100_m}},
      {6, {200_um, 200_um}},
      {4, {200_um, 200_um}},
      {2, {200_um, 200_um}},
      {0, {200_um, 200_um}}};

  LayerPathWidthProvider pathWidthProvider(pathWidths);

  pathSeederCfg.pathWidthProvider.connect<&LayerPathWidthProvider::operator()>(
      &pathWidthProvider);

  E320SourceLinkGridConstructor::Config gridConstructorCfg{
      .bins = std::make_pair(1000, 1), .layers = layerMap};
  gridConstructorCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  auto gridConstructor =
      std::make_shared<E320SourceLinkGridConstructor>(gridConstructorCfg);

  // Estimator of the IP and first hit
  // parameters of the track
  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> refLayers;
  const auto& refVolume = detector->findDetectorVolume("layer0");
  for (const auto* surf : refVolume->surfaces()) {
    refLayers.try_emplace(surf->geometryId(), surf);
    pathSeederCfg.refLayerIds.push_back(surf->geometryId());
  }

  JsonTrackLookupReader::Config lookupReaderCfg;
  lookupReaderCfg.refLayers = refLayers;
  lookupReaderCfg.bins = {1000, 1};

  TrackLookupProvider::Config lookupProviderCfg;
  lookupProviderCfg.lookupPath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_lookups/no_corrector/lookup-prototype.json";
  lookupProviderCfg.trackLookupReader =
      std::make_shared<JsonTrackLookupReader>(lookupReaderCfg);
  TrackLookupProvider lookupProvider(lookupProviderCfg);

  pathSeederCfg.trackEstimator.connect<&TrackLookupProvider::lookup>(
      &lookupProvider);

  // Create the path seeder algorithm
  auto seedingAlgoCfg = PathSeedingAlgorithm::Config();
  seedingAlgoCfg.seeder = std::make_shared<Acts::PathSeeder>(pathSeederCfg);
  seedingAlgoCfg.sourceLinkGridConstructor = gridConstructor;
  seedingAlgoCfg.inputSourceLinks = "Measurements";
  seedingAlgoCfg.outputSeeds = "PathSeeds";
  seedingAlgoCfg.minSeedSize = 1;
  seedingAlgoCfg.maxSeedSize = 1e5;
  sequencer.addAlgorithm(
      std::make_shared<PathSeedingAlgorithm>(seedingAlgoCfg, logLevel));

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
      if (vol->name() == "layer0") {
        cuts.push_back({surf->geometryId(),
                        {{},
                         {std::numeric_limits<double>::max()},
                         {1000u},
                         {std::numeric_limits<double>::max()}}});
      } else {
        /*cuts.push_back({surf->geometryId(), {{}, {15}, {1u}, {15}}});*/
        cuts.push_back({surf->geometryId(),
                        {{},
                         {std::numeric_limits<double>::max()},
                         {1u},
                         {std::numeric_limits<double>::max()}}});
      }
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
  trackFindingCfg.inputSeeds = "PathSeeds";
  trackFindingCfg.outputTrackCandidates = "TrackCandidates";
  trackFindingCfg.outputTrackView = "CandidatesTrackView";
  trackFindingCfg.minCandidateSize = 5;
  trackFindingCfg.maxCandidateSize = 5;
  trackFindingCfg.maxSteps = 1e4;

  auto trackFindingAlgorithm =
      std::make_shared<CKFTrackFindingAlgorithm>(trackFindingCfg, logLevel);
  sequencer.addAlgorithm(trackFindingAlgorithm);

  // --------------------------------------------------------------
  // Alignment

  Acts::GainMatrixUpdater kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;

  // Initialize track fitter options
  Acts::KalmanFitterExtensions<KFTrajectory> alignmentExtensions;
  // Add calibrator
  alignmentExtensions.calibrator
      .connect<&simpleSourceLinkCalibrator<KFTrajectory>>();
  // Add the updater
  alignmentExtensions.updater
      .connect<&Acts::GainMatrixUpdater::operator()<KFTrajectory>>(&kfUpdater);
  // Add the smoother
  alignmentExtensions.smoother
      .connect<&Acts::GainMatrixSmoother::operator()<KFTrajectory>>(
          &kfSmoother);
  // Add the surface accessor
  alignmentExtensions.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  auto alignmentPropOptions = PropagatorOptions(gctx, mctx);

  alignmentPropOptions.maxSteps = 1000;

  auto alignmentKFOptions = Acts::KalmanFitterOptions(
      gctx, mctx, cctx, alignmentExtensions, alignmentPropOptions);

  alignmentKFOptions.referenceSurface = refSurface.get();

  ActsAlignment::AlignedTransformUpdater voidAlignUpdater =
      [&alignCtx](Acts::DetectorElementBase* element,
                  const Acts::GeometryContext& gctx,
                  const Acts::Transform3& transform) {
        std::cout << "\n\n\nUPDATER CALL\n";
        std::cout << "AT: " << element->surface().geometryId() << "\n";
        std::cout << "TRANSLATION: " << transform.translation().transpose()
                  << "\n";
        std::cout << "ROTATION: " << transform.rotation() << "\n";
        alignCtx.alignmentStore->at(element->surface().geometryId()) =
            transform;
        return true;
      };
  AlignmentTransformUpdater transformUpdater;

  AlignmentAlgorithm::Config alignmentCfg{
      .inputTrackCandidates = "TrackCandidates",
      .outputAlignmentParameters = "AlignmentParameters",
      .referenceSurface = refSurface,
      .align = AlignmentAlgorithm::makeAlignmentFunction(detector, field),
      .alignedTransformUpdater = voidAlignUpdater,
      .kfOptions = alignmentKFOptions,
      .chi2ONdfCutOff = 1e-6,
      .maxNumIterations = 3};

  for (auto& det : detector->detectorElements()) {
    const auto& surface = det->surface();
    if (surface.geometryId().sensitive() != 9 &&
        surface.geometryId().sensitive() != 1) {
      alignmentCfg.alignedDetElements.push_back(det.get());
    }
  }

  auto alignmentAlgorithm =
      std::make_shared<AlignmentAlgorithm>(alignmentCfg, logLevel);
  sequencer.addAlgorithm(alignmentAlgorithm);

  // --------------------------------------------------------------
  // Track fitting

  // Initialize track fitter options
  Acts::KalmanFitterExtensions<KFTrajectory> extensions;
  // Add calibrator
  extensions.calibrator.connect<&simpleSourceLinkCalibrator<KFTrajectory>>();
  // Add the updater
  extensions.updater
      .connect<&Acts::GainMatrixUpdater::operator()<KFTrajectory>>(&kfUpdater);
  // Add the smoother
  extensions.smoother
      .connect<&Acts::GainMatrixSmoother::operator()<KFTrajectory>>(
          &kfSmoother);
  // Add the surface accessor
  extensions.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  auto propOptions = PropagatorOptions(gctx, mctx);

  propOptions.maxSteps = 300;

  auto options =
      Acts::KalmanFitterOptions(gctx, mctx, cctx, extensions, propOptions);

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
  TrackFittingAlgorithm::Config fitterCfg{
      .inputTrackCandidates = "TrackCandidates",
      .outputTracks = "Tracks",
      .fitter = fitter,
      .kfOptions = options};

  sequencer.addAlgorithm(
      std::make_shared<TrackFittingAlgorithm>(fitterCfg, logLevel));

  // --------------------------------------------------------------
  // Event writeout

  // Sim cluster writer
  RootSimClusterWriter::Config clusterWriterCfg;

  clusterWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);
  clusterWriterCfg.inputClusters = "Clusters";
  clusterWriterCfg.treeName = "clusters";
  clusterWriterCfg.filePath = "clusters.root";

  sequencer.addWriter(
      std::make_shared<RootSimClusterWriter>(clusterWriterCfg, logLevel));

  // Seed writer
  RootSimSeedWriter::Config seedWriterCfg;

  seedWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);
  seedWriterCfg.inputSeeds = "PathSeeds";
  seedWriterCfg.inputTruthClusters = "Clusters";
  seedWriterCfg.treeName = "seeds";
  seedWriterCfg.filePath = "seeds.root";

  sequencer.addWriter(
      std::make_shared<RootSimSeedWriter>(seedWriterCfg, logLevel));

  // Track candidate writer
  RootSimTrackCandidateWriter::Config trackCandidateWriterCfg;
  trackCandidateWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  trackCandidateWriterCfg.inputTrackCandidates = "CandidatesTrackView";
  trackCandidateWriterCfg.inputTruthClusters = "Clusters";
  trackCandidateWriterCfg.treeName = "track-candidates";
  trackCandidateWriterCfg.filePath = "track-candidates.root";
  trackCandidateWriterCfg.targetTrueTrackSize = 5;

  sequencer.addWriter(std::make_shared<RootSimTrackCandidateWriter>(
      trackCandidateWriterCfg, logLevel));

  //  // Alignment parameters writer
  //  AlignmentParametersWriter::Config alignmentWriterCfg;
  //  alignmentWriterCfg.inputAlignmentResults = "AlignmentParameters";
  //
  //  sequencer.addWriter(std::make_shared<AlignmentParametersWriter>(
  //      alignmentWriterCfg, logLevel));

  // Fitted track writer
  auto trackWriterCfg = RootSimTrackWriter::Config();
  trackWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  trackWriterCfg.inputTracks = "Tracks";
  trackWriterCfg.inputTruthClusters = "Clusters";
  trackWriterCfg.treeName = "fitted-tracks";
  trackWriterCfg.filePath = "fitted-tracks.root";
  trackWriterCfg.targetTrueTrackSize = 5;

  //  sequencer.addWriter(
  //      std::make_shared<RootSimTrackWriter>(trackWriterCfg, logLevel));

  sequencer.run();
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive()) {
        std::cout << s->center(gctx).transpose() << "\n";
      }
    }
  }
  return 0;
}
