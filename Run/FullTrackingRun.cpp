#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Seeding/PathSeeder.hpp"
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

#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/AlignmentParametersProvider.hpp"
#include "TrackingPipeline/Io/E320RootDataReader.hpp"
#include "TrackingPipeline/Io/RootMeasurementWriter.hpp"
#include "TrackingPipeline/Io/RootNoamSplitDataReader.hpp"
#include "TrackingPipeline/Io/RootSeedWriter.hpp"
#include "TrackingPipeline/Io/RootTrackCandidateWriter.hpp"
#include "TrackingPipeline/Io/RootTrackReader.hpp"
#include "TrackingPipeline/Io/RootTrackWriter.hpp"
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

  Acts::Vector3 xCorrectorB(0, 0, -0.026107_T);
  ConstantBoundedField xCorrectorField(xCorrectorB, xCorrectorExtent);

  CompositeMagField::FieldComponents fieldComponents = {
      {quad1Extent, &quad1Field},
      {quad2Extent, &quad2Field},
      {quad3Extent, &quad3Field},
      {dipoleExtent, &dipoleField},
      {xCorrectorExtent, &xCorrectorField}};

  auto field = std::make_shared<CompositeMagField>(fieldComponents);

  auto aStore =
      std::make_shared<std::map<Acts::GeometryIdentifier, Acts::Transform3>>();

  double detectorTilt = 0.0;
  std::map<int, Acts::Vector3> shifts{
      {8, Acts::Vector3(-11.10528_mm + 0.00_um, 0_mm, -(0.7354_mm + 0.00_um))},
      {6, Acts::Vector3(-11.10528_mm - 1.93_um, 0_mm, -(0.7354_mm + 16.5_um))},
      {4, Acts::Vector3(-11.10528_mm + 41.8_um, 0_mm, -(0.7354_mm - 63.3_um))},
      {2, Acts::Vector3(-11.10528_mm + 74.6_um, 0_mm, -(0.7354_mm - 69.5_um))},
      {0, Acts::Vector3(-11.10528_mm - 17.6_um, 0_mm, -(0.7354_mm + 27.5_um))}};
  std::map<int, double> rots{{8, 0_rad},
                             {6, -1.73e-3_rad},
                             {4, 1.46e-3_rad},
                             {2, -1.49e-3_rad},
                             {0, 5.16e-4_rad}};

  // Alignment provider
  auto alignmentProviderCfg = AlignmentParametersProvider::Config();

  alignmentProviderCfg.treeName = "alignment-results";
  alignmentProviderCfg.filePath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_analysis/noam_split/temp/alignment_15/"
      "alignment-results.root";

  // AlignmentParametersProvider alignmentProvider(alignmentProviderCfg);
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      // if (s->geometryId().sensitive() && s->geometryId().sensitive() != 9 &&
      //     s->geometryId().sensitive() != 1) {
      //   Acts::Transform3 nominal = s->transform(Acts::GeometryContext());
      //   const auto [shift, rot] =
      //       alignmentProvider.getAlignedTransform(s->geometryId());
      //   nominal.pretranslate(shifts.at(s->geometryId().sensitive() - 1));
      //   nominal.pretranslate(shift);
      //   nominal.rotate(rot);
      //   aStore->emplace(s->geometryId(), nominal);
      // } else if (s->geometryId().sensitive() == 9 ||
      //            s->geometryId().sensitive() == 1) {
      if (s->geometryId().sensitive()) {
        Acts::Transform3 nominal = s->transform(Acts::GeometryContext());
        nominal.pretranslate(shifts.at(s->geometryId().sensitive() - 1));
        nominal.rotate(Acts::AngleAxis3(
            rots.at(s->geometryId().sensitive() - 1), Acts::Vector3::UnitZ()));
        aStore->emplace(s->geometryId(), nominal);
      }
    }
  }
  AlignmentContext alignCtx(aStore);
  gctx = Acts::GeometryContext{alignCtx};

  Acts::GeometryIdentifier midGeoId;
  midGeoId.setSensitive(5);
  Acts::Vector3 detectorCenter = detector->findSurface(midGeoId)->center(gctx);
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive()) {
        // if (s->geometryId().sensitive() == 9 ||
        //     s->geometryId().sensitive() == 1) {
        Acts::Transform3 nominal = aStore->at(s->geometryId());
        nominal.pretranslate(-detectorCenter);
        nominal.prerotate(Acts::AngleAxis3(detectorTilt, Acts::Vector3::UnitX())
                              .toRotationMatrix());
        nominal.pretranslate(detectorCenter);
        aStore->at(s->geometryId()) = nominal;
      }
    }
  }

  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive()) {
        std::cout << "-----------------------------------\n";
        std::cout << "SURFACE " << s->geometryId() << "\n";
        std::cout << "CENTER " << s->center(gctx).transpose() << " -- "
                  << s->center(Acts::GeometryContext()).transpose() << "\n";
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

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));

  // Add the sim data reader
  E320Io::RootNoamSplitDataReader::Config readerCfg;
  readerCfg.treeName = "MyTree";
  readerCfg.outputSourceLinks = "Measurements";
  readerCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);
  std::string pathToDir =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_dataInRootFormat/"
      "E320Shift_Nov_2024/noam_split/test";

  // Get the paths to the files in the directory
  for (const auto& entry : std::filesystem::directory_iterator(pathToDir)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".root") {
      continue;
    }
    std::string pathToFile = entry.path();
    readerCfg.filePaths.push_back(pathToFile);
  }

  // Add the reader to the sequencer
  // sequencer.addReader(std::make_shared<RootTrackReader>(readerCfg,
  // logLevel));
  sequencer.addReader(
      std::make_shared<E320Io::RootNoamSplitDataReader>(readerCfg, logLevel));

  // --------------------------------------------------------------
  // The path seeding setup
  auto pathSeederCfg = Acts::PathSeeder::Config();

  // Intersection finder
  std::vector<const Acts::Surface*> layers;
  for (const auto& vol : detector->volumes()) {
    for (const auto& surf : vol->surfaces()) {
      if (surf->geometryId().sensitive()) {
        layers.push_back(surf);
      }
    }
  }

  ForwardOrderedIntersectionFinder::Config intersectionFinderCfg;
  intersectionFinderCfg.layers = std::move(layers);
  ForwardOrderedIntersectionFinder intersectionFinder(intersectionFinderCfg);

  pathSeederCfg.intersectionFinder
      .connect<&ForwardOrderedIntersectionFinder::operator()>(
          &intersectionFinder);

  const auto& refVolume = detector->findDetectorVolume("layer0");
  std::vector<Acts::GeometryIdentifier> refGeoIds;
  for (const auto* surf : refVolume->surfaces()) {
    refGeoIds.push_back(surf->geometryId());
  }
  pathSeederCfg.refLayerIds = std::move(refGeoIds);

  // Path width provider
  std::map<int, std::pair<double, double>> pathWidths = {{8, {100_m, 100_m}},
                                                         {6, {150_um, 150_um}},
                                                         {4, {250_um, 250_um}},
                                                         {2, {350_um, 350_um}},
                                                         {0, {500_um, 500_um}}};
  LayerPathWidthProvider pathWidthProvider(pathWidths);

  pathSeederCfg.pathWidthProvider.connect<&LayerPathWidthProvider::operator()>(
      &pathWidthProvider);

  // Grid to bin the source links
  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> layerMap;
  for (const auto& vol : detector->volumes()) {
    for (const auto& surf : vol->surfaces()) {
      if (surf->geometryId().sensitive()) {
        layerMap.insert({surf->geometryId(), surf});
      }
    }
  }

  E320SourceLinkGridConstructor::Config gridConstructorCfg{
      .bins = std::make_pair(250, 125), .layers = layerMap};
  gridConstructorCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  auto gridConstructor =
      std::make_shared<E320SourceLinkGridConstructor>(gridConstructorCfg);

  // Estimator of the IP and first hit
  // parameters of the track
  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> refLayers;
  for (const auto* surf : refVolume->surfaces()) {
    refLayers.try_emplace(surf->geometryId(), surf);
  }

  E320DipoleTrackLookupProvider::Config lookupProviderCfg;
  lookupProviderCfg.dipoleAmplidute = 0.2192;
  lookupProviderCfg.dipolePosition = gOpt.dipoleTranslation[2];
  lookupProviderCfg.dipoleSize = 0.914;

  lookupProviderCfg.correctorAmplidute = -0.026107_T;
  lookupProviderCfg.correctorPosition = gOpt.xCorrectorTranslation[2];
  lookupProviderCfg.correctorSize = 0.23368;

  lookupProviderCfg.layerPosition = gOpt.staveZ.at(8);
  lookupProviderCfg.detectorYZTilt = detectorTilt;
  lookupProviderCfg.referenceSurface = refLayers.begin()->second;
  E320DipoleTrackLookupProvider lookupProvider(lookupProviderCfg);

  pathSeederCfg.trackEstimator.connect<&E320DipoleTrackLookupProvider::lookup>(
      &lookupProvider);

  // Create the path seeder algorithm
  auto seedingAlgoCfg = PathSeedingAlgorithm::Config();
  seedingAlgoCfg.seeder = std::make_shared<Acts::PathSeeder>(pathSeederCfg);
  seedingAlgoCfg.sourceLinkGridConstructor = gridConstructor;
  seedingAlgoCfg.inputSourceLinks = "Measurements";
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
                        {{}, {std::numeric_limits<double>::max()}, {1000u}}});
      } else {
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

  // double refZ = 13060.6_mm + 914_mm / 2;
  double refZ = gOpt.beWindowTranslation[2];
  // double refZ = gOpt.staveZ.at(8) - 0.1_mm;
  Acts::Transform3 transform(Acts::Translation3(Acts::Vector3(0, refZ, 0)) *
                             gOpt.actsToWorldRotation.inverse());

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

  sequencer.addAlgorithm(
      std::make_shared<KFTrackFittingAlgorithm>(fitterCfg, logLevel));

  // --------------------------------------------------------------
  // Event write out

  // Sim cluster writer
  auto clusterWriterCfg = RootMeasurementWriter::Config();

  clusterWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  clusterWriterCfg.inputMeasurements = "Measurements";
  clusterWriterCfg.treeName = "clusters";
  clusterWriterCfg.filePath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/TrackingPipeline_build/"
      "clusters-data.root";

  sequencer.addWriter(
      std::make_shared<RootMeasurementWriter>(clusterWriterCfg, logLevel));

  // Seed writer
  auto seedWriterCfg = RootSeedWriter::Config();

  seedWriterCfg.inputSeeds = "PathSeeds";
  seedWriterCfg.treeName = "seeds";
  seedWriterCfg.filePath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/TrackingPipeline_build/"
      "seeds-data.root";

  seedWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  sequencer.addWriter(
      std::make_shared<RootSeedWriter>(seedWriterCfg, logLevel));

  // Track candidate writer
  auto trackCandidateWriterCfg = RootTrackCandidateWriter::Config();
  trackCandidateWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  trackCandidateWriterCfg.inputTrackCandidates = "CandidatesTrackView";
  trackCandidateWriterCfg.treeName = "track-candidates";
  trackCandidateWriterCfg.filePath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/TrackingPipeline_build/"
      "track-candidates-data.root";

  sequencer.addWriter(std::make_shared<RootTrackCandidateWriter>(
      trackCandidateWriterCfg, logLevel));

  // Fitted track writer
  auto trackWriterCfg = RootTrackWriter::Config();
  trackWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  trackWriterCfg.inputTracks = "Tracks";
  trackWriterCfg.treeName = "fitted-tracks";
  trackWriterCfg.filePath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/TrackingPipeline_build/"
      "fitted-tracks-data.root";

  sequencer.addWriter(
      std::make_shared<RootTrackWriter>(trackWriterCfg, logLevel));

  return sequencer.run();
}
