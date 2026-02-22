#include "TrackingPipeline/Run/SimAlignmentRun.hpp"

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
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <toml.hpp>

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
#include "TrackingPipeline/Io/AlignmentParametersWriter.hpp"
#include "TrackingPipeline/Io/AlignmentParametersProvider.hpp"
#include "TrackingPipeline/Io/RootSimSeedWriter.hpp"
#include "TrackingPipeline/Io/RootSimTrackReader.hpp"
#include "TrackingPipeline/Io/RootSimTrackWriter.hpp"
#include "TrackingPipeline/TrackFitting/KFTrackFittingAlgorithm.hpp"

using namespace Acts::UnitLiterals;
namespace ag = E320Geometry;

// Propagator short-hands
using ActionList = Acts::ActionList<>;
using AbortList  = Acts::AbortList<Acts::EndOfWorldReached>;

using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                    Acts::Experimental::DetectorNavigator>;
using PropagatorOptions =
    typename Propagator::template Options<ActionList, AbortList>;

// KF short-hands
using RecoTrajectory     = KFTrackFittingAlgorithm::Trajectory;
using RecoTrackContainer = KFTrackFittingAlgorithm::TrackContainer;
using KF                 = Acts::KalmanFitter<Propagator, RecoTrajectory>;

std::unique_ptr<const ag::GeometryOptions> ag::GeometryOptions::m_instance =
    nullptr;

namespace TrackingPipeline {

SimAlignmentRunConfig parseSimAlignmentRunConfig(const std::string& path) {
  auto root = toml::parse(path);
  SimAlignmentRunConfig cfg;
  cfg.root = root;

  // [MAIN]
  if (root.contains("MAIN")) {
    const auto& mainTbl = toml::find(root, "MAIN");
    cfg.main.logLevel =
        toml::find_or<std::string>(mainTbl, "log_level", "INFO");
    cfg.main.outputDir =
        toml::find_or<std::string>(mainTbl, "output_dir", std::string{});
  }

  // [SEQUENCER]
  if (root.contains("SEQUENCER")) {
    const auto& seqTbl = toml::find(root, "SEQUENCER");
    cfg.sequencer.numThreads =
        toml::find_or<std::size_t>(seqTbl, "num_threads", 1u);
    cfg.sequencer.skip =
        toml::find_or<std::size_t>(seqTbl, "skip", 0u);
    cfg.sequencer.trackFpes =
        toml::find_or<bool>(seqTbl, "track_fpes", false);
  }

  // [TRACK_READER]
  {
    const auto& rTbl = toml::find(root, "TRACK_READER");
    cfg.trackReader.inputFiles =
        toml::find<std::vector<std::string>>(rTbl, "input_files");
    cfg.trackReader.treeName =
        toml::find_or<std::string>(rTbl, "tree_name", "fitted-tracks");
    cfg.trackReader.outputMeasurements =
        toml::find_or<std::string>(rTbl, "output_measurements", "SimMeasurements");
    cfg.trackReader.outputSimClusters =
        toml::find_or<std::string>(rTbl, "output_sim_clusters", "SimClusters");
    cfg.trackReader.outputSeedsGuess =
        toml::find_or<std::string>(rTbl, "output_seeds_guess", "SeedsGuess");
    cfg.trackReader.outputSeedsEst =
        toml::find_or<std::string>(rTbl, "output_seeds_est", "SeedsEst");
    cfg.trackReader.minChi2 =
        toml::find_or<double>(rTbl, "min_chi2", 0.0);
    cfg.trackReader.maxChi2 =
        toml::find_or<double>(rTbl, "max_chi2", 1e4);
    cfg.trackReader.mergeIntoOneEvent =
        toml::find_or<bool>(rTbl, "merge_into_one_event", true);
  }

  // [ALIGNMENT_PROVIDER]
  if (root.contains("ALIGNMENT_PROVIDER")) {
  const auto& aTbl = toml::find(root, "ALIGNMENT_PROVIDER");
  cfg.alignmentProvider.useFile =
      toml::find_or<bool>(aTbl, "use_file", false);
  cfg.alignmentProvider.filePath =
      toml::find_or<std::string>(aTbl, "file_path", std::string{});
  cfg.alignmentProvider.treeName =
      toml::find_or<std::string>(aTbl, "tree_name", "alignment-parameters");
  }

  // [ALIGNMENT_KF]
  if (root.contains("ALIGNMENT_KF")) {
  const auto& kTbl = toml::find(root, "ALIGNMENT_KF");
  cfg.alignmentKF.maxSteps =
      toml::find_or<std::size_t>(kTbl, "max_steps", 1000u);

  double loc0_mm       = toml::find_or<double>(kTbl, "loc0_mm", 100.0);
  double loc1_mm       = toml::find_or<double>(kTbl, "loc1_mm", 100.0);
  double phi_rad       = toml::find_or<double>(kTbl, "phi_rad", 10.0);
  double theta_rad     = toml::find_or<double>(kTbl, "theta_rad", 10.0);
  double qoverp_GeVinv = toml::find_or<double>(kTbl, "q_over_p_GeVinv", 100.0);
  double time_ns       = toml::find_or<double>(kTbl, "time_ns", 25.0);

  cfg.alignmentKF.originStdDevLoc0   = loc0_mm   * 1_mm;
  cfg.alignmentKF.originStdDevLoc1   = loc1_mm   * 1_mm;
  cfg.alignmentKF.originStdDevPhi    = phi_rad   * 1_rad;
  cfg.alignmentKF.originStdDevTheta  = theta_rad * 1_rad;
  cfg.alignmentKF.originStdDevQOverP = qoverp_GeVinv / 1_GeV;
  cfg.alignmentKF.originStdDevTime   = time_ns   * 1_ns;
  } 

  // [ALIGNMENT_ALGO]
  if (root.contains("ALIGNMENT_ALGO")) {
    const auto& aTbl = toml::find(root, "ALIGNMENT_ALGO");
    cfg.alignmentAlgo.inputTrackCandidates =
        toml::find_or<std::string>(aTbl, "input_track_candidates", "SeedsEst");
    cfg.alignmentAlgo.outputAlignmentParameters =
        toml::find_or<std::string>(aTbl, "output_alignment_parameters", "AlignmentParameters");
    cfg.alignmentAlgo.chi2ONdfCutOff =
        toml::find_or<double>(aTbl, "chi2_on_ndf_cutoff", 1e-16);
    cfg.alignmentAlgo.deltaChi2ONdfCutOff0 =
        toml::find_or<double>(aTbl, "delta_chi2_on_ndf_cutoff0", 10.0);
    cfg.alignmentAlgo.deltaChi2ONdfCutOff1 =
        toml::find_or<double>(aTbl, "delta_chi2_on_ndf_cutoff1", 1e-5);
    cfg.alignmentAlgo.maxNumIterations =
        toml::find_or<std::size_t>(aTbl, "max_num_iterations", 200u);
    cfg.alignmentAlgo.alignmentMode =
        toml::find_or<std::string>(aTbl, "alignment_mode", "global");
    cfg.alignmentAlgo.propagationDirection =
        toml::find_or<std::string>(aTbl, "propagation_direction", "backward");
  }

  // [FITTING]
  if (root.contains("FITTING")) {
    const auto& fTbl = toml::find(root, "FITTING");
    cfg.fitting.maxSteps =
        toml::find_or<std::size_t>(fTbl, "max_steps", 100000u);
    cfg.fitting.resolvePassive =
        toml::find_or<bool>(fTbl, "resolve_passive", false);
    cfg.fitting.resolveMaterial =
        toml::find_or<bool>(fTbl, "resolve_material", true);
    cfg.fitting.resolveSensitive =
        toml::find_or<bool>(fTbl, "resolve_sensitive", true);
    cfg.fitting.inputCandidates =
        toml::find_or<std::string>(fTbl, "input_candidates", "SeedsGuess");
    cfg.fitting.outputTracks =
        toml::find_or<std::string>(fTbl, "output_tracks", "Tracks");
  }

  // [SEED_WRITER]
  if (root.contains("SEED_WRITER")) {
    const auto& sTbl = toml::find(root, "SEED_WRITER");
    cfg.seedWriter.enable =
        toml::find_or<bool>(sTbl, "enable", true);
    cfg.seedWriter.inputSeeds =
        toml::find_or<std::string>(sTbl, "inputSeeds", "SeedsEst");
    cfg.seedWriter.inputSimClusters =
        toml::find_or<std::string>(sTbl, "inputSimClusters", "SimClusters");
    cfg.seedWriter.treeName =
        toml::find_or<std::string>(sTbl, "treeName", "seeds");
    cfg.seedWriter.outputFile =
        toml::find_or<std::string>(sTbl, "output_file", "seeds.root");
  }

  // [TRACK_WRITER]
  if (root.contains("TRACK_WRITER")) {
    const auto& tTbl = toml::find(root, "TRACK_WRITER");
    cfg.trackWriter.enable =
        toml::find_or<bool>(tTbl, "enable", true);
    cfg.trackWriter.inputTracks =
        toml::find_or<std::string>(tTbl, "inputTracks", "Tracks");
    cfg.trackWriter.inputSimClusters =
        toml::find_or<std::string>(tTbl, "inputSimClusters", "SimClusters");
    cfg.trackWriter.treeName =
        toml::find_or<std::string>(tTbl, "treeName", "fitted-tracks");
    cfg.trackWriter.outputFile =
        toml::find_or<std::string>(tTbl, "output_file", "fitted-tracks.root");
  }

  // [ALIGNMENT_WRITER]
  if (root.contains("ALIGNMENT_WRITER")) {
    const auto& aWTbl = toml::find(root, "ALIGNMENT_WRITER");
    cfg.alignmentWriter.enable =
        toml::find_or<bool>(aWTbl, "enable", true);
    cfg.alignmentWriter.inputAlignmentParameters =
        toml::find_or<std::string>(aWTbl, "inputAlignmentParameters", "AlignmentParameters");
    cfg.alignmentWriter.treeName =
        toml::find_or<std::string>(aWTbl, "treeName", "alignment-parameters");
    cfg.alignmentWriter.outputFile =
        toml::find_or<std::string>(aWTbl, "output_file", "alignment-parameters.root");
  }

  // [CONSTRAINTS] (optional)
  if (root.contains("CONSTRAINTS")) {
    const auto& cTbl = toml::find(root, "CONSTRAINTS");
    cfg.constraints.minSensitiveIdForConstraints =
        toml::find_or<int>(cTbl, "min_sensitive_id_for_constraints", 40);
    cfg.constraints.minSensitiveIdForAlign =
        toml::find_or<int>(cTbl, "min_sensitive_id_for_align", 10);
    cfg.constraints.maxSensitiveIdForAlign =
        toml::find_or<int>(cTbl, "max_sensitive_id_for_align", 40);
  }

  return cfg;
}

Acts::Logging::Level getLogLevel(const std::string& levelStr) {
  if (levelStr == "VERBOSE") return Acts::Logging::VERBOSE;
  if (levelStr == "DEBUG")   return Acts::Logging::DEBUG;
  if (levelStr == "INFO")    return Acts::Logging::INFO;
  if (levelStr == "WARNING") return Acts::Logging::WARNING;
  if (levelStr == "ERROR")   return Acts::Logging::ERROR;
  if (levelStr == "FATAL")   return Acts::Logging::FATAL;
  return Acts::Logging::INFO;
}

int runSimAlignment(const std::string& configPath) {
  SimAlignmentRunConfig cfg;
  try {
    cfg = parseSimAlignmentRunConfig(configPath);
  } catch (const std::exception& e) {
    std::cerr << "Error parsing SimAlignmentRun config: " << e.what() << "\n";
    return 1;
  }

  Acts::Logging::Level logLevel = getLogLevel(cfg.main.logLevel);
  const auto& goInst            = *ag::GeometryOptions::instance();

  // Contexts
  Acts::GeometryContext      gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext   cctx;

  // Detector
  auto detector = E320Geometry::buildDetector(gctx);

  // Alignment store: from file or built
  std::shared_ptr<AlignmentContext::AlignmentStore> aStore;
  if (cfg.alignmentProvider.useFile &&
      !cfg.alignmentProvider.filePath.empty()) {
    AlignmentParametersProvider::Config alignmentProviderCfg;
    alignmentProviderCfg.filePath = cfg.alignmentProvider.filePath;
    alignmentProviderCfg.treeName = cfg.alignmentProvider.treeName;
    AlignmentParametersProvider alignmentProvider(alignmentProviderCfg);
    aStore = alignmentProvider.getAlignmentStore();
  } else {
    aStore = detail::makeAlignmentStore(gctx, detector.get());
  }

  AlignmentContext alignCtx(aStore);
  Acts::GeometryContext testCtx{alignCtx};
  gctx = Acts::GeometryContext{alignCtx};

  // Magnetic field
  auto field = E320Geometry::buildMagField(gctx);

  // --------------------------------------------------------------
  // Event reading
  SimpleSourceLink::SurfaceAccessor   simpleSurfaceAccessor{detector.get()};
  ExtendedSourceLink::SurfaceAccessor extendedSurfaceAccessor{detector.get()};

  MixedSourceLinkSurfaceAccessor surfaceAccessor;
  surfaceAccessor.connect<&SimpleSourceLink::SurfaceAccessor::operator(),
                          SimpleSourceLink>(&simpleSurfaceAccessor);
  surfaceAccessor.connect<&ExtendedSourceLink::SurfaceAccessor::operator(),
                          ExtendedSourceLink>(&extendedSurfaceAccessor);

  // Sequencer
  Sequencer::Config seqCfg;
  seqCfg.numThreads = cfg.sequencer.numThreads;
  seqCfg.skip       = cfg.sequencer.skip;
  seqCfg.trackFpes  = cfg.sequencer.trackFpes;
  Sequencer sequencer(seqCfg);

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));

  // Track reader
  RootSimTrackReader::Config readerCfg;
  readerCfg.treeName          = cfg.trackReader.treeName;
  readerCfg.outputMeasurements = cfg.trackReader.outputMeasurements;
  readerCfg.outputSimClusters  = cfg.trackReader.outputSimClusters;
  readerCfg.outputSeedsGuess   = cfg.trackReader.outputSeedsGuess;
  readerCfg.outputSeedsEst     = cfg.trackReader.outputSeedsEst;
  readerCfg.minChi2            = cfg.trackReader.minChi2;
  readerCfg.maxChi2            = cfg.trackReader.maxChi2;
  readerCfg.mergeIntoOneEvent  = cfg.trackReader.mergeIntoOneEvent;

  for (const auto& path : cfg.trackReader.inputFiles) {
    readerCfg.filePaths.push_back(path);
  }

  sequencer.addReader(
      std::make_shared<RootSimTrackReader>(readerCfg, logLevel));

  // Reference surface
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

  // Alignment constraints
  std::vector<Acts::GeometryIdentifier> constraintsSurfaceIds;
  for (auto& det : detector->detectorElements()) {
    const auto& surface = det->surface();
    const auto& sid     = surface.geometryId().sensitive();
    if (sid && sid >= cfg.constraints.minSensitiveIdForConstraints) {
      constraintsSurfaceIds.push_back(surface.geometryId());
    }
  }

  std::vector<Acts::SourceLink> alignmentConstraints;
  for (const auto& gid : constraintsSurfaceIds) {
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
        ExtendedSourceLink(loc, glob, cov, gid, 0, 0));
  }

  // Calibrators
  MixedSourceLinkCalibrator<RecoTrajectory> mixedSourceLinkCalibrator;
  mixedSourceLinkCalibrator.connect<
      &extendedSourceLinkCalibrator<RecoTrajectory>, ExtendedSourceLink>();
  mixedSourceLinkCalibrator
      .connect<&simpleSourceLinkCalibrator<RecoTrajectory>, SimpleSourceLink>();

  // --------------------------------------------------------------
  // Alignment KF setup
  Acts::GainMatrixUpdater  kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;

  Acts::KalmanFitterExtensions<RecoTrajectory> alignmentExtensions;
  alignmentExtensions.calibrator
      .connect<&MixedSourceLinkCalibrator<RecoTrajectory>::operator()>(
          &mixedSourceLinkCalibrator);
  alignmentExtensions.updater
      .connect<&Acts::GainMatrixUpdater::operator()<RecoTrajectory>>(
          &kfUpdater);
  alignmentExtensions.smoother
      .connect<&Acts::GainMatrixSmoother::operator()<RecoTrajectory>>(
          &kfSmoother);
  alignmentExtensions.surfaceAccessor
      .connect<&MixedSourceLinkSurfaceAccessor::operator()>(&surfaceAccessor);

  auto alignmentPropOptions = PropagatorOptions(gctx, mctx);
  alignmentPropOptions.maxSteps = cfg.alignmentKF.maxSteps;

  auto alignmentKFOptions = Acts::KalmanFitterOptions(
      gctx, mctx, cctx, alignmentExtensions, alignmentPropOptions);
  alignmentKFOptions.referenceSurface = refSurface.get();

  // Initial track state covariance matrix
  Acts::BoundVector trackOriginStdDevPrior;
  trackOriginStdDevPrior[Acts::eBoundLoc0]   = cfg.alignmentKF.originStdDevLoc0;
  trackOriginStdDevPrior[Acts::eBoundLoc1]   = cfg.alignmentKF.originStdDevLoc1;
  trackOriginStdDevPrior[Acts::eBoundTime]   = cfg.alignmentKF.originStdDevTime;
  trackOriginStdDevPrior[Acts::eBoundPhi]    = cfg.alignmentKF.originStdDevPhi;
  trackOriginStdDevPrior[Acts::eBoundTheta]  = cfg.alignmentKF.originStdDevTheta;
  trackOriginStdDevPrior[Acts::eBoundQOverP] = cfg.alignmentKF.originStdDevQOverP;
  Acts::BoundMatrix trackOriginCov =
      trackOriginStdDevPrior.cwiseProduct(trackOriginStdDevPrior).asDiagonal();

  AlignmentAlgorithm::Config alignmentCfg{
      .inputTrackCandidates      = cfg.alignmentAlgo.inputTrackCandidates,
      .outputAlignmentParameters = cfg.alignmentAlgo.outputAlignmentParameters,
      .align                     = AlignmentAlgorithm::makeAlignmentFunction(detector, field),
      .alignedTransformUpdater   = detail::makeGlobalAlignmentUpdater(alignCtx),
      .kfOptions                 = alignmentKFOptions,
      .chi2ONdfCutOff            = cfg.alignmentAlgo.chi2ONdfCutOff,
      .deltaChi2ONdfCutOff       = {cfg.alignmentAlgo.deltaChi2ONdfCutOff0, cfg.alignmentAlgo.deltaChi2ONdfCutOff1},
      .maxNumIterations          = cfg.alignmentAlgo.maxNumIterations,
      .alignmentMask             = (ActsAlignment::AlignmentMask::Center1 |
                                    ActsAlignment::AlignmentMask::Center2 |
                                    ActsAlignment::AlignmentMask::Rotation2),
      .alignmentMode             = ActsAlignment::AlignmentMode::global,
      .originCov                 = trackOriginCov,
      .constraints               = alignmentConstraints,
      .propDirection             = AlignmentAlgorithm::PropagationDirection::backward};

  if (cfg.alignmentAlgo.alignmentMode == "global")
    alignmentCfg.alignmentMode = ActsAlignment::AlignmentMode::global;
  else if (cfg.alignmentAlgo.alignmentMode == "local")
    alignmentCfg.alignmentMode = ActsAlignment::AlignmentMode::local;

  if (cfg.alignmentAlgo.propagationDirection == "forward")
    alignmentCfg.propDirection = AlignmentAlgorithm::PropagationDirection::forward;
  else
    alignmentCfg.propDirection = AlignmentAlgorithm::PropagationDirection::backward;

  for (auto& det : detector->detectorElements()) {
    const auto& surface = det->surface();
    const auto sid      = surface.geometryId().sensitive();
    if (sid &&
        sid >= cfg.constraints.minSensitiveIdForAlign &&
        sid <  cfg.constraints.maxSensitiveIdForAlign) {
      alignmentCfg.alignedDetElements.push_back(det.get());
    }
  }

  auto alignmentAlgorithm =
      std::make_shared<AlignmentAlgorithm>(alignmentCfg, logLevel);
  sequencer.addAlgorithm(alignmentAlgorithm);

  // --------------------------------------------------------------
  // Track fitting KF setup
  Acts::KalmanFitterExtensions<RecoTrajectory> extensions;
  extensions.calibrator.connect<&simpleSourceLinkCalibrator<RecoTrajectory>>();
  extensions.updater
      .connect<&Acts::GainMatrixUpdater::operator()<RecoTrajectory>>(
          &kfUpdater);
  extensions.smoother
      .connect<&Acts::GainMatrixSmoother::operator()<RecoTrajectory>>(
          &kfSmoother);
  extensions.surfaceAccessor
      .connect<&MixedSourceLinkSurfaceAccessor::operator()>(&surfaceAccessor);

  auto propOptions = PropagatorOptions(gctx, mctx);
  propOptions.maxSteps = cfg.fitting.maxSteps;

  auto kfOptions = Acts::KalmanFitterOptions(
      gctx, mctx, cctx, extensions, propOptions, refSurface.get());

  Acts::Experimental::DetectorNavigator::Config navCfg;
  navCfg.detector        = detector.get();
  navCfg.resolvePassive  = cfg.fitting.resolvePassive;
  navCfg.resolveMaterial = cfg.fitting.resolveMaterial;
  navCfg.resolveSensitive = cfg.fitting.resolveSensitive;
  Acts::Experimental::DetectorNavigator kfNavigator(
      navCfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));

  Acts::EigenStepper<> kfStepper(std::move(field));
  auto kfPropagator =
      Propagator(std::move(kfStepper), std::move(kfNavigator),
                 Acts::getDefaultLogger("Propagator", logLevel));

  const auto fitter = KF(
      kfPropagator, Acts::getDefaultLogger("DetectorKalmanFilter", logLevel));

  KFTrackFittingAlgorithm::Config fitterCfg{
      .inputTrackCandidates = cfg.fitting.inputCandidates,
      .outputTracks         = cfg.fitting.outputTracks,
      .fitter               = fitter,
      .kfOptions            = kfOptions};

  sequencer.addAlgorithm(
      std::make_shared<KFTrackFittingAlgorithm>(fitterCfg, logLevel));

  // --------------------------------------------------------------
  // Writers
  std::filesystem::path outDir(cfg.main.outputDir);
  try {
    std::filesystem::create_directories(outDir);
  } catch (const std::exception& e) {
    std::cerr << "Error creating output directory '" << outDir.string()
              << "': " << e.what() << "\n";
    return 1;
  }

  // Seed writer
  if (cfg.seedWriter.enable) {
    RootSimSeedWriter::Config seedWriterCfg;
    seedWriterCfg.inputSeeds        = cfg.seedWriter.inputSeeds;
    seedWriterCfg.inputTruthClusters = cfg.seedWriter.inputSimClusters;
    seedWriterCfg.treeName          = cfg.seedWriter.treeName;
    seedWriterCfg.filePath          =
        (outDir / cfg.seedWriter.outputFile).string();

    sequencer.addWriter(
        std::make_shared<RootSimSeedWriter>(seedWriterCfg, logLevel));
  }

  // Track writer
  if (cfg.trackWriter.enable) {
    RootSimTrackWriter::Config trackWriterCfg;
    trackWriterCfg.surfaceAccessor
        .connect<&MixedSourceLinkSurfaceAccessor::operator()>(&surfaceAccessor);
    trackWriterCfg.referenceSurface = refSurface.get();
    trackWriterCfg.inputTracks      = cfg.trackWriter.inputTracks;
    trackWriterCfg.inputSimClusters = cfg.trackWriter.inputSimClusters;
    trackWriterCfg.treeName         = cfg.trackWriter.treeName;
    trackWriterCfg.filePath         =
        (outDir / cfg.trackWriter.outputFile).string();

    sequencer.addWriter(
        std::make_shared<RootSimTrackWriter>(trackWriterCfg, logLevel));
  }

  // Alignment parameters writer
  if (cfg.alignmentWriter.enable) {
    AlignmentParametersWriter::Config alignmentWriterCfg;
    alignmentWriterCfg.treeName              = cfg.alignmentWriter.treeName;
    alignmentWriterCfg.inputAlignmentResults = cfg.alignmentWriter.inputAlignmentParameters;
    alignmentWriterCfg.filePath              =
        (outDir / cfg.alignmentWriter.outputFile).string();

    sequencer.addWriter(
        std::make_shared<AlignmentParametersWriter>(alignmentWriterCfg, logLevel));
  }

  sequencer.run();
  return 0;
}

} // namespace TrackingPipeline

int main(int argc, char* argv[]) {
  std::string confPath = (argc > 1)
                           ? argv[1]
                           : "SimAlignmentRun.conf";
  return TrackingPipeline::runSimAlignment(confPath);
}
