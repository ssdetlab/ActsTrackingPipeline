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
#include <unordered_map>

#include <Eigen/src/Core/ArithmeticSequence.h>
#include <nlohmann/json.hpp>
#include <unistd.h>

#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"
#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Alignment/detail/AlignmentUtils.hpp"
#include "TrackingPipeline/Geometry/ApollonGeometry.hpp"
#include "TrackingPipeline/Geometry/ApollonGeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/AlignmentParametersProvider.hpp"
#include "TrackingPipeline/Io/AlignmentParametersWriter.hpp"
#include "TrackingPipeline/Io/RootSeedWriter.hpp"
#include "TrackingPipeline/Io/RootTrackReader.hpp"
#include "TrackingPipeline/Io/RootTrackWriter.hpp"
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
      std::cout << surf->polyhedronRepresentation(gctx, 1000).extent() << "\n";
      if (surf->geometryId().sensitive()) {
        surfaceMap[surf->geometryId()] = surf;
      }
    }
  }

//   AlignmentParametersProvider::Config alignmentProviderCfg1;
//   alignmentProviderCfg1.filePath =
//       "/home/romanurmanov/work/Apollon/tracking/out_data/Apollon_cosmic_data/"
//       "alignment/stave_1/aligned_1/"
//       "alignment-parameters.root";
//   alignmentProviderCfg1.treeName = "alignment-parameters";
//   AlignmentParametersProvider alignmentProvider1(alignmentProviderCfg1);
//   auto aStore1 = alignmentProvider1.getAlignmentStore();

//   AlignmentParametersProvider::Config alignmentProviderCfg2;
//   alignmentProviderCfg2.filePath =
//       "/home/romanurmanov/work/Apollon/tracking/out_data/Apollon_cosmic_data/"
//       "alignment/stave_0/aligned/"
//       "alignment-parameters.root";
//   alignmentProviderCfg2.treeName = "alignment-parameters";
//   AlignmentParametersProvider alignmentProvider2(alignmentProviderCfg2);
//   auto aStore2 = alignmentProvider2.getAlignmentStore();

//   auto aStore = std::make_shared<AlignmentContext::AlignmentStore>();
//   for (const auto& entry : *aStore1) {
//     aStore->insert(entry);
//   }
//   for (const auto& entry : *aStore2) {
//     aStore->insert(entry);
//   }

//   AlignmentContext alignCtx(aStore);
//   Acts::GeometryContext testCtx{alignCtx};
//   for (auto& v : detector->volumes()) {
//     for (auto& s : v->surfaces()) {
//       if (s->geometryId().sensitive()) {
//         std::cout << "-----------------------------------\n";
//         std::cout << "SURFACE " << s->geometryId() << "\n";
//         std::cout << "CENTER " << s->center(testCtx).transpose() << " -- "
//                   << s->center(Acts::GeometryContext()).transpose() << "\n";
//         std::cout << "NORMAL "
//                   << s->normal(testCtx, s->center(testCtx),
//                                Acts::Vector3::UnitY())
//                          .transpose()
//                   << " -- "
//                   << s->normal(testCtx, s->center(Acts::GeometryContext()),
//                                Acts::Vector3::UnitY())
//                          .transpose()
//                   << "\n";
//         std::cout << "ROTATION \n"
//                   << s->transform(testCtx).rotation() << " -- \n"
//                   << "\n"
//                   << s->transform(Acts::GeometryContext()).rotation() << "\n";

//         std::cout << "EXTENT "
//                   << s->polyhedronRepresentation(testCtx, 1000).extent()
//                   << "\n -- \n"
//                   << s->polyhedronRepresentation(Acts::GeometryContext(), 1000)
//                          .extent()
//                   << "\n";
//       }
//     }
//   }
//   gctx = Acts::GeometryContext{alignCtx};

  double longTransStd = 0_um;
  double shortTransStd = 0_um;
  std::size_t longIdx = detail::binningValueToIndex(goInst.longBinValue);
  std::size_t shortIdx = detail::binningValueToIndex(goInst.shortBinValue);
  auto aStore = detail::makeAlignmentStore(
      detector.get(), longIdx, longTransStd, shortIdx, shortTransStd);
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
                  << s->polyhedronRepresentation(Acts::GeometryContext(),
                  1000)
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
  seqCfg.numThreads = 1;
  seqCfg.skip = 0;
  seqCfg.trackFpes = false;
  Sequencer sequencer(seqCfg);

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));

  // Add the sim data reader
  RootTrackReader::Config readerCfg;
  readerCfg.treeName = "fitted-tracks";
  readerCfg.outputMeasurements = "SimMeasurements";
  readerCfg.outputSeeds = "Seeds";
  readerCfg.minChi2 = 250;
  readerCfg.maxChi2 = 650;
  readerCfg.mergeIntoOneEvent = true;

  std::string pathToDir =
      "/Users/nathalyn/Desktop/Masters/ACTS-data/output/filtered";

  // Get the paths to the files in the directory
  for (const auto& entry : std::filesystem::directory_iterator(pathToDir)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".root") {
      continue;
    }
    std::string pathToFile = entry.path();
    readerCfg.filePaths.push_back(pathToFile);
  }

  // Add the reader to the sequencer
  sequencer.addReader(std::make_shared<RootTrackReader>(readerCfg, logLevel));

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

  Acts::Transform3 transform = Acts::Transform3::Identity();
  transform.translate(Acts::Vector3(0_mm, 0, 0));
  transform.rotate(refSurfToWorldRotationX);
  transform.rotate(refSurfToWorldRotationY);
  transform.rotate(refSurfToWorldRotationZ);

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
  Acts::KalmanFitterExtensions<RecoTrajectory> alignmentExtensions;
  // Add calibrator
  alignmentExtensions.calibrator
      .connect<&simpleSourceLinkCalibrator<RecoTrajectory>>();
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
        alignCtx.alignmentStore()[element->surface().geometryId()] = transform;
        return true;
      };
  AlignmentTransformUpdater transformUpdater;

  AlignmentAlgorithm::Config alignmentCfg{
      .inputTrackCandidates = "Seeds",
      .outputAlignmentParameters = "AlignmentParameters",
      .align = AlignmentAlgorithm::makeAlignmentFunction(detector, field),
      .alignedTransformUpdater = voidAlignUpdater,
      .kfOptions = alignmentKFOptions,
      .chi2ONdfCutOff = 1e-3,
      .deltaChi2ONdfCutOff = {10, 1e-5},
      .maxNumIterations = 200,
      .alignmentMask = (ActsAlignment::AlignmentMask::Center1 |
                        ActsAlignment::AlignmentMask::Center2 |
                        ActsAlignment::AlignmentMask::Rotation2),
      .alignmentMode = ActsAlignment::AlignmentMode::local};

  for (auto& det : detector->detectorElements()) {
    const auto& surface = det->surface();
    const auto& geoId = surface.geometryId().sensitive();
    if (geoId &&
        // surface.geometryId().sensitive() >= 22) {
        surface.geometryId().sensitive() < 22 &&
        surface.geometryId().sensitive() != 10) {
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
  KFTrackFittingAlgorithm::Config fitterCfg{.inputTrackCandidates = "Seeds",
                                            .outputTracks = "Tracks",
                                            .fitter = fitter,
                                            .kfOptions = options};

  sequencer.addAlgorithm(
      std::make_shared<KFTrackFittingAlgorithm>(fitterCfg, logLevel));

  // --------------------------------------------------------------
  // Event write out

  // Seed writer
  RootSeedWriter::Config seedWriterCfg;
  seedWriterCfg.inputSeeds = "Seeds";
  seedWriterCfg.inputTruthClusters = "SimClusters";
  seedWriterCfg.treeName = "seeds";
  seedWriterCfg.filePath =
      "/Users/nathalyn/Desktop/Masters/ACTS-data/output/seeds.root";

  sequencer.addWriter(
      std::make_shared<RootSeedWriter>(seedWriterCfg, logLevel));

  // Fitted track writer
  RootTrackWriter::Config trackWriterCfg;
  trackWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);
  trackWriterCfg.referenceSurface = refSurface.get();
  trackWriterCfg.inputTracks = "Tracks";
  trackWriterCfg.treeName = "fitted-tracks";
  trackWriterCfg.filePath =
      "/Users/nathalyn/Desktop/Masters/ACTS-data/output/"
      "fitted-tracks.root";

  sequencer.addWriter(
      std::make_shared<RootTrackWriter>(trackWriterCfg, logLevel));

  // Alignment parameters writer
  AlignmentParametersWriter::Config alignmentWriterCfg;
  alignmentWriterCfg.treeName = "alignment-parameters";
  alignmentWriterCfg.inputAlignmentResults = "AlignmentParameters";
  alignmentWriterCfg.filePath =
      "/Users/nathalyn/Desktop/Masters/ACTS-data/output/"
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
                  << s->normal(gctx, s->center(Acts::GeometryContext()),
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