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

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include <Eigen/src/Core/ArithmeticSequence.h>
#include <nlohmann/json.hpp>
#include <unistd.h>

#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"
#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Alignment/detail/AlignmentStoreBuilders.hpp"
#include "TrackingPipeline/Alignment/detail/AlignmentStoreUpdaterBuilders.hpp"
#include "TrackingPipeline/EventData/ExtendedSourceLink.hpp"
#include "TrackingPipeline/EventData/MixedSourceLinkCalibrator.hpp"
#include "TrackingPipeline/EventData/MixedSourceLinkSurfaceAccessor.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/AlignmentParametersProvider.hpp"
#include "TrackingPipeline/Io/AlignmentParametersWriter.hpp"
#include "TrackingPipeline/Io/RootSimSeedWriter.hpp"
#include "TrackingPipeline/Io/RootSimTrackReader.hpp"
#include "TrackingPipeline/Io/RootSimTrackWriter.hpp"
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

namespace ag = E320Geometry;

std::unique_ptr<const ag::GeometryOptions> ag::GeometryOptions::m_instance =
    nullptr;

int main() {
  const auto& goInst = *ag::GeometryOptions::instance();

  // Set the log level
  Acts::Logging::Level logLevel = Acts::Logging::INFO;

  // Dummy context and options
  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext cctx;

  // --------------------------------------------------------------
  // Detector setup

  auto detector = E320Geometry::buildDetector(gctx);

  // AlignmentParametersProvider::Config alignmentProviderCfg;
  // alignmentProviderCfg.filePath =
  //     "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/sim/"
  //     "alignment/global/aligned/"
  //     "alignment-parameters.root";
  // alignmentProviderCfg.treeName = "alignment-parameters";
  // AlignmentParametersProvider alignmentProvider(alignmentProviderCfg);
  // auto aStore = alignmentProvider.getAlignmentStore();

  // AlignmentContext alignCtx(aStore);
  // Acts::GeometryContext testCtx{alignCtx};
  // for (auto& v : detector->volumes()) {
  //   for (auto& s : v->surfaces()) {
  //     if (s->geometryId().sensitive()) {
  //       std::cout << "-----------------------------------\n";
  //       std::cout << "SURFACE " << s->geometryId() << "\n";
  //       std::cout << "CENTER " << s->center(testCtx).transpose() << " -- "
  //                 << s->center(Acts::GeometryContext()).transpose() << "\n";
  //       std::cout << "NORMAL "
  //                 << s->normal(testCtx, s->center(testCtx),
  //                              Acts::Vector3::UnitY())
  //                        .transpose()
  //                 << " -- "
  //                 << s->normal(testCtx, s->center(Acts::GeometryContext()),
  //                              Acts::Vector3::UnitY())
  //                        .transpose()
  //                 << "\n";
  //       std::cout << "ROTATION \n"
  //                 << s->transform(testCtx).rotation() << " -- \n"
  //                 << "\n"
  //                 << s->transform(Acts::GeometryContext()).rotation() << "\n";

  //       std::cout << "EXTENT "
  //                 << s->polyhedronRepresentation(testCtx, 1000).extent()
  //                 << "\n -- \n"
  //                 << s->polyhedronRepresentation(Acts::GeometryContext(),
  //                 1000)
  //                        .extent()
  //                 << "\n";
  //     }
  //   }
  // }
  // gctx = Acts::GeometryContext{alignCtx};

  auto aStore = detail::makeAlignmentStore(gctx, detector.get());
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

  auto field = E320Geometry::buildMagField(gctx);

  // --------------------------------------------------------------
  // Event reading
  SimpleSourceLink::SurfaceAccessor simpleSurfaceAccessor{detector.get()};
  ExtendedSourceLink::SurfaceAccessor extendedSurfaceAccessor{detector.get()};

  MixedSourceLinkSurfaceAccessor surfaceAccessor;
  surfaceAccessor.connect<&SimpleSourceLink::SurfaceAccessor::operator(),
                          SimpleSourceLink>(&simpleSurfaceAccessor);
  surfaceAccessor.connect<&ExtendedSourceLink::SurfaceAccessor::operator(),
                          ExtendedSourceLink>(&extendedSurfaceAccessor);

  // Setup the sequencer
  Sequencer::Config seqCfg;
  // seqCfg.events = 1e1;
  seqCfg.numThreads = 1;
  seqCfg.skip = 0;
  seqCfg.trackFpes = false;
  Sequencer sequencer(seqCfg);

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));

  // Add the sim data reader
  RootSimTrackReader::Config readerCfg;
  readerCfg.treeName = "fitted-tracks";
  readerCfg.outputMeasurements = "SimMeasurements";
  readerCfg.outputSimClusters = "SimClusters";
  readerCfg.outputSeedsGuess = "SeedsGuess";
  readerCfg.outputSeedsEst = "SeedsEst";
  readerCfg.minChi2 = 0;
  readerCfg.maxChi2 = 1e4;
  readerCfg.mergeIntoOneEvent = true;

  std::string pathToDir =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/sim/"
      "alignment/global/misaligned/fitted_tracks";

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
      std::make_shared<RootSimTrackReader>(readerCfg, logLevel));

  // --------------------------------------------------------------
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

  Acts::Transform3 refSurfaceTransform = Acts::Transform3::Identity();
  refSurfaceTransform.translate(Acts::Vector3(0, 0, 0));
  refSurfaceTransform.rotate(refSurfToWorldRotationX);
  refSurfaceTransform.rotate(refSurfToWorldRotationY);
  refSurfaceTransform.rotate(refSurfToWorldRotationZ);

  auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      refSurfaceTransform,
      std::make_shared<Acts::RectangleBounds>(halfX, halfY));

  Acts::GeometryIdentifier geoId;
  geoId.setExtra(1);
  refSurface->assignGeometryId(std::move(geoId));

  // --------------------------------------------------------------
  // Alignment

  // Initialize constraints
  std::vector<Acts::GeometryIdentifier> constraintsSurfaceIds;
  for (auto& det : detector->detectorElements()) {
    const auto& surface = det->surface();
    const auto& geoId = surface.geometryId();
    if (geoId.sensitive() && geoId.sensitive() >= 40) {
      constraintsSurfaceIds.push_back(geoId);
    }
  }

  std::vector<Acts::SourceLink> alignmentConstraints;
  for (const auto& geoId : constraintsSurfaceIds) {
    Acts::ActsVector<7> glob =
        Acts::ActsVector<ExtendedSourceLink::globalSubspaceSize>::Zero();
    Acts::ActsVector<ExtendedSourceLink::localSubspaceSize> loc =
        Acts::ActsVector<ExtendedSourceLink::localSubspaceSize>::Zero();
    loc(3) = M_PI_2;
    loc(4) = 1.0 / 2.5_GeV;
    Acts::ActsVector<ExtendedSourceLink::localSubspaceSize> stdDev = {
        10_mm, 10_mm, 10_rad, 10_rad, 1 / 0.1_GeV};
    Acts::ActsSquareMatrix<ExtendedSourceLink::localSubspaceSize> cov =
        stdDev.cwiseProduct(stdDev).asDiagonal();
    alignmentConstraints.emplace_back(
        ExtendedSourceLink(loc, glob, cov, geoId, 0, 0));
  }

  // Initialize calibrators
  MixedSourceLinkCalibrator<RecoTrajectory> mixedSourceLinkCalibrator;
  mixedSourceLinkCalibrator.connect<
      &extendedSourceLinkCalibrator<RecoTrajectory>, ExtendedSourceLink>();
  mixedSourceLinkCalibrator
      .connect<&simpleSourceLinkCalibrator<RecoTrajectory>, SimpleSourceLink>();

  // Initialize track fitter options
  Acts::GainMatrixUpdater kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;

  Acts::KalmanFitterExtensions<RecoTrajectory> alignmentExtensions;
  // Add calibrator
  alignmentExtensions.calibrator
      .connect<&MixedSourceLinkCalibrator<RecoTrajectory>::operator()>(
          &mixedSourceLinkCalibrator);
  // Add the updater
  alignmentExtensions.updater
      .connect<&Acts::GainMatrixUpdater::operator()<RecoTrajectory>>(
          &kfUpdater);
  // Add the smoother
  alignmentExtensions.smoother
      .connect<&Acts::GainMatrixSmoother::operator()<RecoTrajectory>>(
          &kfSmoother);
  // Add the surface accessor
  alignmentExtensions.surfaceAccessor
      .connect<&MixedSourceLinkSurfaceAccessor::operator()>(&surfaceAccessor);

  auto alignmentPropOptions = PropagatorOptions(gctx, mctx);

  alignmentPropOptions.maxSteps = 1000;

  auto alignmentKFOptions = Acts::KalmanFitterOptions(
      gctx, mctx, cctx, alignmentExtensions, alignmentPropOptions);

  alignmentKFOptions.referenceSurface = refSurface.get();

  // Initial track state covariance matrix
  Acts::BoundVector trackOriginStdDevPrior;
  trackOriginStdDevPrior[Acts::eBoundLoc0] = 100_mm;
  trackOriginStdDevPrior[Acts::eBoundLoc1] = 100_mm;
  trackOriginStdDevPrior[Acts::eBoundTime] = 25_ns;
  trackOriginStdDevPrior[Acts::eBoundPhi] = 10_rad;
  trackOriginStdDevPrior[Acts::eBoundTheta] = 10_rad;
  trackOriginStdDevPrior[Acts::eBoundQOverP] = 1 / 0.01_GeV;
  Acts::BoundMatrix trackOriginCov =
      trackOriginStdDevPrior.cwiseProduct(trackOriginStdDevPrior).asDiagonal();

  AlignmentAlgorithm::Config alignmentCfg{
      .inputTrackCandidates = "SeedsEst",
      .outputAlignmentParameters = "AlignmentParameters",
      .align = AlignmentAlgorithm::makeAlignmentFunction(detector, field),
      .alignedTransformUpdater = detail::makeGlobalAlignmentUpdater(alignCtx),
      .kfOptions = alignmentKFOptions,
      .chi2ONdfCutOff = 1e-16,
      .deltaChi2ONdfCutOff = {10, 1e-5},
      .maxNumIterations = 200,
      .alignmentMask = (ActsAlignment::AlignmentMask::Center1 |
                        ActsAlignment::AlignmentMask::Center2 |
                        ActsAlignment::AlignmentMask::Rotation2),
      .alignmentMode = ActsAlignment::AlignmentMode::global,
      .originCov = trackOriginCov,
      .constraints = alignmentConstraints,
      .propDirection = AlignmentAlgorithm::PropagationDirection::backward};

  for (auto& det : detector->detectorElements()) {
    const auto& surface = det->surface();
    const auto& geoId = surface.geometryId().sensitive();
    if (geoId && surface.geometryId().sensitive() >= 10 &&
        surface.geometryId().sensitive() < 40) {
      alignmentCfg.alignedDetElements.push_back(det.get());
    }
  }

  auto alignmentAlgorithm =
      std::make_shared<AlignmentAlgorithm>(alignmentCfg, logLevel);
  sequencer.addAlgorithm(alignmentAlgorithm);

  // --------------------------------------------------------------
  // Track fitting

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
      .connect<&MixedSourceLinkSurfaceAccessor::operator()>(&surfaceAccessor);

  auto propOptions = PropagatorOptions(gctx, mctx);

  propOptions.maxSteps = 1e5;

  auto kfOptions = Acts::KalmanFitterOptions(gctx, mctx, cctx, extensions,
                                             propOptions, refSurface.get());

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
      .inputTrackCandidates = "SeedsGuess",
      .outputTracks = "Tracks",
      .fitter = fitter,
      .kfOptions = kfOptions};

  sequencer.addAlgorithm(
      std::make_shared<KFTrackFittingAlgorithm>(fitterCfg, logLevel));

  // --------------------------------------------------------------
  // Event write out

  // Seed writer
  RootSimSeedWriter::Config seedWriterCfg;
  seedWriterCfg.inputSeeds = "SeedsEst";
  seedWriterCfg.inputTruthClusters = "SimClusters";
  seedWriterCfg.treeName = "seeds";
  seedWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/sim/"
      "seeds.root";

  sequencer.addWriter(
      std::make_shared<RootSimSeedWriter>(seedWriterCfg, logLevel));

  // Fitted track writer
  RootSimTrackWriter::Config trackWriterCfg;
  trackWriterCfg.surfaceAccessor
      .connect<&MixedSourceLinkSurfaceAccessor::operator()>(&surfaceAccessor);
  trackWriterCfg.referenceSurface = refSurface.get();
  trackWriterCfg.inputTracks = "Tracks";
  trackWriterCfg.inputSimClusters = "SimClusters";
  trackWriterCfg.treeName = "fitted-tracks";
  trackWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/sim/"
      "fitted-tracks.root";

  sequencer.addWriter(
      std::make_shared<RootSimTrackWriter>(trackWriterCfg, logLevel));

  // Alignment parameters writer
  AlignmentParametersWriter::Config alignmentWriterCfg;
  alignmentWriterCfg.treeName = "alignment-parameters";
  alignmentWriterCfg.inputAlignmentResults = "AlignmentParameters";
  alignmentWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_analysis/sim/"
      "alignment-parameters.root";

  sequencer.addWriter(std::make_shared<AlignmentParametersWriter>(
      alignmentWriterCfg, logLevel));

  sequencer.run();

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

  return 0;
}
