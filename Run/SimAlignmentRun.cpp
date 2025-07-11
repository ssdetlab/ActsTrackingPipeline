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

#include <unistd.h>

#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"
#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/RootSimClusterReader.hpp"
#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"
#include "TrackingPipeline/Io/RootSimSeedWriter.hpp"
#include "TrackingPipeline/Io/RootSimTrackCandidateWriter.hpp"
#include "TrackingPipeline/Io/RootSimTrackReader.hpp"
#include "TrackingPipeline/Io/RootSimTrackWriter.hpp"
#include "TrackingPipeline/MagneticField/CompositeMagField.hpp"
#include "TrackingPipeline/MagneticField/ConstantBoundedField.hpp"
#include "TrackingPipeline/MagneticField/DipoleMagField.hpp"
#include "TrackingPipeline/MagneticField/IdealQuadrupoleMagField.hpp"
#include "TrackingPipeline/TrackFinding/CKFTrackFindingAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/DipoleTrackLookupProvider.hpp"
#include "TrackingPipeline/TrackFinding/E320SourceLinkGridConstructor.hpp"
#include "TrackingPipeline/TrackFinding/ForwardOrderedIntersectionFinder.hpp"
#include "TrackingPipeline/TrackFinding/LayerPathWidthProvider.hpp"
#include "TrackingPipeline/TrackFinding/PathSeedingAlgorithm.hpp"
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

int main() {
  // Set the log level
  Acts::Logging::Level logLevel = Acts::Logging::DEBUG;

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
  std::map<int, Acts::Vector3> shifts{
      {8, Acts::Vector3(-9.7_mm, -3.5_mm, 0_mm)},
      {6, Acts::Vector3(-9.7_mm, -3.5_mm, 0_mm)},
      {4, Acts::Vector3(-9.7_mm, -3.5_mm, 0_mm)},
      {2, Acts::Vector3(-9.7_mm, -3.5_mm, 0_mm)},
      {0, Acts::Vector3(-9.7_mm, -3.5_mm, 0_mm)}};
  Acts::RotationMatrix3 mat8 =
      Acts::AngleAxis3(0, Acts::Vector3::UnitZ()).toRotationMatrix();
  Acts::RotationMatrix3 mat6 =
      Acts::AngleAxis3(0, Acts::Vector3::UnitZ()).toRotationMatrix();
  Acts::RotationMatrix3 mat4 =
      Acts::AngleAxis3(0, Acts::Vector3::UnitZ()).toRotationMatrix();
  Acts::RotationMatrix3 mat2 =
      Acts::AngleAxis3(0, Acts::Vector3::UnitZ()).toRotationMatrix();
  Acts::RotationMatrix3 mat0 =
      Acts::AngleAxis3(0, Acts::Vector3::UnitZ()).toRotationMatrix();

  std::map<int, Acts::RotationMatrix3> rots{
      {8, mat8}, {6, mat6}, {4, mat4}, {2, mat2}, {0, mat0}};
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive()) {
        // Surface is in origin, normal along z, no rotation in xy plane
        Acts::Transform3 nominal = Acts::Transform3::Identity();

        // Global detector rotation
        nominal.rotate(gOpt.actsToWorldRotation.inverse());

        // Global detector translation
        nominal.translate(
            Acts::Vector3(gOpt.chipX, gOpt.chipY,
                          gOpt.staveZ.at(s->geometryId().sensitive() - 1)));

        // Apply relative translations of the rotated surfaces
        nominal.translate(shifts.at(s->geometryId().sensitive() - 1));

        // Rotate surface in the origin around global origin
        nominal.rotate(rots.at(s->geometryId().sensitive() - 1));

        // Account for G4 rotation
        nominal.rotate(Acts::AngleAxis3(-M_PI_2, Acts::Vector3::UnitZ())
                           .toRotationMatrix());

        aStore->emplace(s->geometryId(), nominal);
      }
    }
  }

  AlignmentContext alignCtx(aStore);
  Acts::GeometryContext testCtx{alignCtx};
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive()) {
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

  IdealQuadrupoleMagField quad1Field(
      gOpt.quadrupolesParams[0],
      gOpt.actsToWorldRotation.inverse() * gOpt.quad1Translation,
      gOpt.actsToWorldRotation);
  IdealQuadrupoleMagField quad2Field(
      gOpt.quadrupolesParams[1],
      gOpt.actsToWorldRotation.inverse() * gOpt.quad2Translation,
      gOpt.actsToWorldRotation);
  IdealQuadrupoleMagField quad3Field(
      gOpt.quadrupolesParams[2],
      gOpt.actsToWorldRotation.inverse() * gOpt.quad3Translation,
      gOpt.actsToWorldRotation);

  double dipoleB = 0.2192_T;
  DipoleMagField dipoleField(
      gOpt.dipoleParams, dipoleB, gOpt.actsToWorldRotation,
      gOpt.actsToWorldRotation.inverse() * gOpt.dipoleTranslation);

  // TODO: Add the real field value
  Acts::Vector3 xCorrectorB(0, 0, -0.026107_T);
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
  seqCfg.events = 1e1;
  seqCfg.numThreads = 1;
  seqCfg.trackFpes = false;
  Sequencer sequencer(seqCfg);

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));

  // --------------------------------------------------------------
  // Add sim clusters reader
  // RootSimTrackReader::Config readerCfg;
  // readerCfg.filePaths = {
  //     "/home/romanurmanov/lab/LUXE/acts_tracking/TrackingPipeline_build/"
  //     "fitted-tracks.root"};
  // readerCfg.treeName = "fitted-tracks";
  // readerCfg.outputMeasurements = "SimMeasurements";
  // readerCfg.outputClusters = "SimClusters";
  // readerCfg.outputSeeds = "PreFittedTrackCandidates";
  // readerCfg.batch = false;
  // readerCfg.batchSize = 100;
  // readerCfg.covAnnealingFactor = 1e0;

  // sequencer.addReader(
  //     std::make_shared<RootSimTrackReader>(readerCfg, logLevel));

  RootSimClusterReader::Config readerCfg;
  readerCfg.filePaths = {
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_analysis/sim/misaligned_local_global/clusters_sim/"
      "clusters-sim.root"};
  readerCfg.treeName = "clusters";
  readerCfg.outputMeasurements = "SimMeasurements";
  readerCfg.outputSimClusters = "SimClusters";

  sequencer.addReader(
      std::make_shared<RootSimClusterReader>(readerCfg, logLevel));

  // --------------------------------------------------------------
  // Reference surface for sampling the track at the IP
  //
  double halfX = std::numeric_limits<double>::max();
  double halfY = std::numeric_limits<double>::max();

  double refZ = gOpt.beWindowTranslation[2] - 2_mm;
  Acts::Transform3 transform(Acts::Translation3(Acts::Vector3(0, refZ, 0)) *
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

  // Nominal
  // std::map<int, std::pair<double, double>> pathWidths = {{8, {100_m, 100_m}},
  //                                                        {6, {250_um,
  //                                                        250_um}}, {4,
  //                                                        {350_um, 350_um}},
  //                                                        {2, {450_um,
  //                                                        450_um}}, {0,
  //                                                        {600_um, 600_um}}};

  // Aligned
  std::map<int, std::pair<double, double>> pathWidths = {{8, {100_m, 100_m}},
                                                         {6, {100_um, 100_um}},
                                                         {4, {100_um, 100_um}},
                                                         {2, {100_um, 100_um}},
                                                         {0, {100_um, 100_um}}};

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

  E320DipoleTrackLookupProvider::Config lookupProviderCfg;
  lookupProviderCfg.dipoleAmplidute = 0.2192;
  lookupProviderCfg.dipolePosition = gOpt.dipoleTranslation[2];
  lookupProviderCfg.dipoleSize = 0.914;

  lookupProviderCfg.correctorAmplidute = -0.026107_T;
  lookupProviderCfg.correctorPosition = gOpt.xCorrectorTranslation[2];
  lookupProviderCfg.correctorSize = 0.23622;

  lookupProviderCfg.layerPosition = gOpt.staveZ.at(8);
  lookupProviderCfg.referenceSurface = refLayers.begin()->second;
  E320DipoleTrackLookupProvider lookupProvider(lookupProviderCfg);

  pathSeederCfg.trackEstimator.connect<&E320DipoleTrackLookupProvider::lookup>(
      &lookupProvider);

  // Create the path seeder algorithm
  auto seedingAlgoCfg = PathSeedingAlgorithm::Config();
  seedingAlgoCfg.seeder = std::make_shared<Acts::PathSeeder>(pathSeederCfg);
  seedingAlgoCfg.sourceLinkGridConstructor = gridConstructor;
  seedingAlgoCfg.inputSourceLinks = "SimMeasurements";
  seedingAlgoCfg.outputSeeds = "PathSeeds";
  seedingAlgoCfg.minSeedSize = 5;
  seedingAlgoCfg.maxSeedSize = 1e5;
  seedingAlgoCfg.minLayers = 5;
  seedingAlgoCfg.maxLayers = 5;
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
                         {100000u},
                         {std::numeric_limits<double>::max()}}});
      } else {
        cuts.push_back({surf->geometryId(),
                        {{},
                         {std::numeric_limits<double>::max()},
                         {2u},
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

  //  // --------------------------------------------------------------
  //  // Alignment
  //
  //  Acts::GainMatrixUpdater kfUpdater;
  //  Acts::GainMatrixSmoother kfSmoother;
  //
  //  // Initialize track fitter options
  //  Acts::KalmanFitterExtensions<KFTrajectory> alignmentExtensions;
  //  // Add calibrator
  //  alignmentExtensions.calibrator
  //      .connect<&simpleSourceLinkCalibrator<KFTrajectory>>();
  //  // Add the updater
  //  alignmentExtensions.updater
  //      .connect<&Acts::GainMatrixUpdater::operator()<KFTrajectory>>(&kfUpdater);
  //  // Add the smoother
  //  alignmentExtensions.smoother
  //      .connect<&Acts::GainMatrixSmoother::operator()<KFTrajectory>>(
  //          &kfSmoother);
  //  // Add the surface accessor
  //  alignmentExtensions.surfaceAccessor
  //      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
  //          &surfaceAccessor);
  //
  //  auto alignmentPropOptions = PropagatorOptions(gctx, mctx);
  //
  //  alignmentPropOptions.maxSteps = 1000;
  //
  //  auto alignmentKFOptions = Acts::KalmanFitterOptions(
  //      gctx, mctx, cctx, alignmentExtensions, alignmentPropOptions);
  //
  //  alignmentKFOptions.referenceSurface = refSurface.get();
  //
  //  ActsAlignment::AlignedTransformUpdater voidAlignUpdater =
  //      [&alignCtx](Acts::DetectorElementBase* element,
  //                  const Acts::GeometryContext& gctx,
  //                  const Acts::Transform3& transform) {
  //        std::cout << "\n\n\nUPDATER CALL\n";
  //        std::cout << "AT: " << element->surface().geometryId() << "\n";
  //        std::cout << "TRANSLATION: " << transform.translation().transpose()
  //                  << "\n";
  //        std::cout << "ROTATION: " << transform.rotation() << "\n";
  //        alignCtx.alignmentStore->at(element->surface().geometryId()) =
  //            transform;
  //        return true;
  //      };
  //  AlignmentTransformUpdater transformUpdater;
  //
  //  AlignmentAlgorithm::Config alignmentCfg{
  //      .inputTrackCandidates = "PreFittedTrackCandidates",
  //      .outputAlignmentParameters = "AlignmentParameters",
  //      .referenceSurface =
  //          detector->findDetectorVolume("beWindow")->surfaces().at(0),
  //      .align = AlignmentAlgorithm::makeAlignmentFunction(detector, field),
  //      .alignedTransformUpdater = voidAlignUpdater,
  //      .kfOptions = alignmentKFOptions,
  //      .chi2ONdfCutOff = 1e-6,
  //      .maxNumIterations = 10};
  //
  //  for (auto& det : detector->detectorElements()) {
  //    const auto& surface = det->surface();
  //    if (surface.geometryId().sensitive() != 9) {
  //      /*if (surface.geometryId().sensitive()) {*/
  //      alignmentCfg.alignedDetElements.push_back(det.get());
  //    }
  //  }
  //
  //  auto alignmentAlgorithm =
  //      std::make_shared<AlignmentAlgorithm>(alignmentCfg, logLevel);
  //  /*sequencer.addAlgorithm(alignmentAlgorithm);*/

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

  propOptions.maxSteps = 1000;

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
  KFTrackFittingAlgorithm::Config fitterCfg{
      .inputTrackCandidates = "TrackCandidates",
      .outputTracks = "Tracks",
      .fitter = fitter,
      .kfOptions = options};

  sequencer.addAlgorithm(
      std::make_shared<KFTrackFittingAlgorithm>(fitterCfg, logLevel));

  // --------------------------------------------------------------
  // Event writeout

  // Sim cluster writer
  RootSimClusterWriter::Config clusterWriterCfg;

  clusterWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);
  clusterWriterCfg.inputClusters = "SimClusters";
  clusterWriterCfg.treeName = "clusters";
  clusterWriterCfg.filePath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_analysis/sim/clusters.root";

  sequencer.addWriter(
      std::make_shared<RootSimClusterWriter>(clusterWriterCfg, logLevel));

  // Seed writer
  RootSimSeedWriter::Config seedWriterCfg;

  seedWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);
  seedWriterCfg.inputSeeds = "PathSeeds";
  seedWriterCfg.inputTruthClusters = "SimClusters";
  seedWriterCfg.treeName = "seeds";
  seedWriterCfg.filePath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_analysis/sim/seeds.root";

  sequencer.addWriter(
      std::make_shared<RootSimSeedWriter>(seedWriterCfg, logLevel));

  // Track candidate writer
  RootSimTrackCandidateWriter::Config trackCandidateWriterCfg;
  trackCandidateWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  trackCandidateWriterCfg.inputTrackCandidates = "CandidatesTrackView";
  trackCandidateWriterCfg.inputTruthClusters = "SimClusters";
  trackCandidateWriterCfg.treeName = "track-candidates";
  trackCandidateWriterCfg.filePath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_analysis/sim/track-candidates.root";

  sequencer.addWriter(std::make_shared<RootSimTrackCandidateWriter>(
      trackCandidateWriterCfg, logLevel));

  // Fitted track writer
  auto trackWriterCfg = RootSimTrackWriter::Config();
  trackWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  trackWriterCfg.inputTracks = "Tracks";
  trackWriterCfg.inputTruthClusters = "SimClusters";
  trackWriterCfg.treeName = "fitted-tracks";
  trackWriterCfg.filePath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_analysis/sim/fitted-tracks.root";

  sequencer.addWriter(
      std::make_shared<RootSimTrackWriter>(trackWriterCfg, logLevel));

  return sequencer.run();
}
