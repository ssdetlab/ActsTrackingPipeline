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
#include "TrackingPipeline/Io/AlignmentParametersProvider.hpp"
#include "TrackingPipeline/Io/AlignmentParametersWriter.hpp"
#include "TrackingPipeline/Io/RootTrackReader.hpp"
#include "TrackingPipeline/Io/RootTrackWriter.hpp"
#include "TrackingPipeline/MagneticField/CompositeMagField.hpp"
#include "TrackingPipeline/MagneticField/ConstantBoundedField.hpp"
#include "TrackingPipeline/MagneticField/DipoleMagField.hpp"
#include "TrackingPipeline/MagneticField/QuadrupoleMagField.hpp"
#include "TrackingPipeline/TrackFitting/KFTrackFittingAlgorithm.hpp"

// Propagator short-hands
using ActionList = Acts::ActionList<>;
using AbortList = Acts::AbortList<Acts::EndOfWorldReached>;

using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                    Acts::Experimental::DetectorNavigator>;
using PropagatorOptions =
    typename Propagator::template Options<ActionList, AbortList>;

// KF short-hands
using KFTrajectory = KFTrackFittingAlgorithm::Trajectory;
using KFTrackContainer = KFTrackFittingAlgorithm::TrackContainer;
using KF = Acts::KalmanFitter<Propagator, KFTrajectory>;

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

  double detectorTilt = 0.0;
  std::map<int, Acts::Vector3> shifts{
      {8, Acts::Vector3(-9.7_mm, 0, -3.5_mm)},
      {6, Acts::Vector3(-9.7_mm, 0, -3.5_mm)},
      {4, Acts::Vector3(-9.7_mm, 0, -3.5_mm)},
      {2, Acts::Vector3(-9.7_mm, 0, -3.5_mm)},
      {0, Acts::Vector3(-9.7_mm, 0, -3.5_mm)}};
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive()) {
        Acts::Transform3 nominal = s->transform(gctx);
        nominal.pretranslate(shifts.at(s->geometryId().sensitive() - 1));
        std::cout << nominal.translation().transpose() << "\n";
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
        Acts::Transform3 nominal = aStore->at(s->geometryId());
        nominal.pretranslate(-detectorCenter);
        nominal.prerotate(Acts::AngleAxis3(detectorTilt, Acts::Vector3::UnitX())
                              .toRotationMatrix());
        nominal.pretranslate(detectorCenter);
        aStore->at(s->geometryId()) = nominal;
      }
    }
  }

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
  seqCfg.events = 1e0;
  seqCfg.numThreads = 1;
  seqCfg.trackFpes = false;
  Sequencer sequencer(seqCfg);

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));

  // --------------------------------------------------------------
  // Add sim clusters reader
  RootTrackReader::Config readerCfg;
  readerCfg.treeName = "fitted-tracks";
  readerCfg.outputMeasurements = "Measurements";
  readerCfg.outputSeeds = "PreFittedTrackCandidates";
  readerCfg.batch = false;
  readerCfg.stack = true;
  readerCfg.batchSize = 1e3;
  readerCfg.covAnnealingFactor = 1e0;
  readerCfg.filePaths = {
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_analysis/data/noam_split/global_yz_tilt_scan/"
      "yz_tilt_test/"
      "initial_full_tracking_run/fitted-tracks-data.root"};

  // // Get the paths to the files in the directory
  // for (const auto& entry : std::filesystem::directory_iterator(pathToDir)) {
  //   if (!entry.is_regular_file() || entry.path().extension() != ".root") {
  //     continue;
  //   }
  //   std::string pathToFile = entry.path();
  //   readerCfg.filePaths.push_back(pathToFile);
  // }
  sequencer.addReader(std::make_shared<RootTrackReader>(readerCfg, logLevel));

  // --------------------------------------------------------------
  // Reference surface for sampling the track at the IP
  double halfX = std::numeric_limits<double>::max();
  double halfY = std::numeric_limits<double>::max();

  double refZ = gOpt.beWindowTranslation[2];
  Acts::Transform3 transform(Acts::Translation3(Acts::Vector3(0, refZ, 0)) *
                             gOpt.actsToWorldRotation.inverse());

  auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      transform, std::make_shared<Acts::RectangleBounds>(halfX, halfY));
  Acts::GeometryIdentifier geoId;
  geoId.setExtra(1);
  refSurface->assignGeometryId(std::move(geoId));

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
        std::cout << element->surface()
                         .polyhedronRepresentation(
                             Acts::GeometryContext(alignCtx), 1000)
                         .extent();
        alignCtx.alignmentStore->at(element->surface().geometryId()) =
            transform;
        return true;
      };
  AlignmentTransformUpdater transformUpdater;

  AlignmentAlgorithm::Config alignmentCfg{
      .inputTrackCandidates = "PreFittedTrackCandidates",
      .outputAlignmentParameters = "AlignmentParameters",
      .referenceSurface =
          detector->findDetectorVolume("beWindow")->surfaces().at(0),
      .align = AlignmentAlgorithm::makeAlignmentFunction(detector, field),
      .alignedTransformUpdater = voidAlignUpdater,
      .kfOptions = alignmentKFOptions,
      .chi2ONdfCutOff = 1e-6,
      .maxNumIterations = 10};

  for (auto& det : detector->detectorElements()) {
    const auto& surface = det->surface();
    if (surface.geometryId().sensitive() != 9) {
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

  propOptions.maxSteps = 1e5;

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
      .inputTrackCandidates = "PreFittedTrackCandidates",
      .outputTracks = "Tracks",
      .fitter = fitter,
      .kfOptions = options};

  sequencer.addAlgorithm(
      std::make_shared<KFTrackFittingAlgorithm>(fitterCfg, logLevel));

  // --------------------------------------------------------------
  // Event writeout

  // Fitted track writer
  auto trackWriterCfg = RootTrackWriter::Config();
  trackWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  trackWriterCfg.inputTracks = "Tracks";
  trackWriterCfg.treeName = "fitted-tracks";
  trackWriterCfg.filePath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_analysis/data/noam_split/fitted-tracks-aligned.root";

  sequencer.addWriter(
      std::make_shared<RootTrackWriter>(trackWriterCfg, logLevel));

  // Alignment writer
  auto alignmentWriterCfg = AlignmentParametersWriter::Config();

  alignmentWriterCfg.inputAlignmentResults = "AlignmentParameters";
  alignmentWriterCfg.treeName = "alignment-results";
  alignmentWriterCfg.filePath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_analysis/data/noam_split/alignment-results.root";

  sequencer.addWriter(std::make_shared<AlignmentParametersWriter>(
      alignmentWriterCfg, logLevel));

  sequencer.run();

  // Alignment provider
  auto alignmentProviderCfg = AlignmentParametersProvider::Config();

  alignmentProviderCfg.treeName = "alignment-results";
  alignmentProviderCfg.filePath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_analysis/data/noam_split/alignment-results.root";

  AlignmentParametersProvider alignmentProvider(alignmentProviderCfg);
  auto bStore =
      std::make_shared<std::map<Acts::GeometryIdentifier, Acts::Transform3>>();
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive() && s->geometryId().sensitive() != 9) {
        Acts::Transform3 nominal = s->transform(Acts::GeometryContext());
        const auto [shift, rot] =
            alignmentProvider.getAlignedTransform(s->geometryId());
        nominal.pretranslate(shift);
        nominal.rotate(rot);
        bStore->emplace(s->geometryId(), nominal);
      } else if (s->geometryId().sensitive() == 9) {
        Acts::Transform3 nominal = s->transform(Acts::GeometryContext());
        nominal.pretranslate(shifts.at(s->geometryId().sensitive() - 1));
        bStore->emplace(s->geometryId(), nominal);
      }
    }
  }
  AlignmentContext blignCtx(bStore);
  Acts::GeometryContext testCtx{blignCtx};

  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive() == 9) {
        Acts::Transform3 nominal = bStore->at(s->geometryId());
        nominal.pretranslate(-detectorCenter);
        nominal.prerotate(Acts::AngleAxis3(detectorTilt, Acts::Vector3::UnitX())
                              .toRotationMatrix());
        nominal.pretranslate(detectorCenter);
        bStore->at(s->geometryId()) = nominal;
      }
    }
  }
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive()) {
        std::cout << "-----------------------------------\n";
        std::cout << "SURFACE " << s->geometryId() << "\n";
        std::cout << "CENTER " << s->center(testCtx).transpose() << " -- "
                  << s->center(gctx).transpose() << "\n";
        std::cout << "NORMAL "
                  << s->normal(testCtx, s->center(testCtx),
                               Acts::Vector3::UnitY())
                         .transpose()
                  << " -- "
                  << s->normal(testCtx, s->center(gctx), Acts::Vector3::UnitY())
                         .transpose()
                  << "\n";
        std::cout << "ROTATION \n"
                  << s->transform(testCtx).rotation() << " -- \n"
                  << "\n"
                  << s->transform(gctx).rotation() << "\n";
        std::cout << "EXTENT "
                  << s->polyhedronRepresentation(testCtx, 1000).extent()
                  << "\n -- \n"
                  << s->polyhedronRepresentation(gctx, 1000).extent() << "\n";
      }
    }
  }
  return 0;
}
