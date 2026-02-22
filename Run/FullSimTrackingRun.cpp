#include "TrackingPipeline/Run/FullSimTrackingRun.hpp"

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

#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <toml.hpp>
#include <vector>

#include <unistd.h>

#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Alignment/detail/AlignmentStoreBuilders.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/AlignmentParametersProvider.hpp"
#include "TrackingPipeline/Io/RootSimClusterReader.hpp"
#include "TrackingPipeline/Io/RootSimSeedWriter.hpp"
#include "TrackingPipeline/Io/RootSimTrackWriter.hpp"
#include "TrackingPipeline/TrackFinding/E320SeedingAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"
#include "TrackingPipeline/TrackFitting/KFTrackFittingAlgorithm.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Io/RootMeasurementWriter.hpp"

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

FullSimTrackingRunConfig parseFullSimTrackingRunConfig(const std::string& path) {
  auto root = toml::parse(path);
  FullSimTrackingRunConfig cfg;
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

  // [READER]
  {
    const auto& inTbl = toml::find(root, "READER");
    cfg.reader.clustersDir =
        toml::find<std::string>(inTbl, "clusters_dir");
    cfg.reader.treeName =
        toml::find_or<std::string>(inTbl, "tree_name", "clusters");
    cfg.reader.minGeoId =
        toml::find_or<int>(inTbl, "min_geo_id", 10);
    cfg.reader.maxGeoId =
        toml::find_or<int>(inTbl, "max_geo_id", 18);
    cfg.reader.surfaceLocalToGlobal =
        toml::find_or<bool>(inTbl, "surface_local_to_global", true);
    cfg.reader.outputSourceLinks =
        toml::find_or<std::string>(inTbl, "outputSourceLinks", "Measurements");
    cfg.reader.outputSimClusters =
        toml::find_or<std::string>(inTbl, "outputSimClusters", "SimClusters");
  }

  // [ALIGNMENT]
  if (root.contains("ALIGNMENT")) {
    const auto& aTbl = toml::find(root, "ALIGNMENT");
    cfg.alignment.useFile =
        toml::find_or<bool>(aTbl, "use_file", false);
    cfg.alignment.filePath =
        toml::find_or<std::string>(aTbl, "file_path", std::string{});
    cfg.alignment.treeName =
        toml::find_or<std::string>(aTbl, "tree_name", "alignment-parameters");
  }

  // [SEEDING]
  if (root.contains("SEEDING")) {
    const auto& sTbl = toml::find(root, "SEEDING");
    cfg.seeding.inputSourceLinks =
        toml::find_or<std::string>(sTbl, "input_source_links", "Measurements");
    cfg.seeding.outputSeeds =
        toml::find_or<std::string>(sTbl, "output_seeds", "Seeds");
    cfg.seeding.minLayers =
        toml::find_or<int>(sTbl, "minLayers", 5);
    cfg.seeding.maxLayers =
        toml::find_or<int>(sTbl, "maxLayers", 5);
    cfg.seeding.beamlineTilt =
        toml::find_or<double>(sTbl, "beamlineTilt", 0.0);

    double loc0_mm        = toml::find_or<double>(sTbl, "loc0_mm", 10.0);
    double loc1_mm        = toml::find_or<double>(sTbl, "loc1_mm", 10.0);
    double phi_deg        = toml::find_or<double>(sTbl, "phi_deg", 10.0);
    double theta_deg      = toml::find_or<double>(sTbl, "theta_deg", 10.0);
    double qoverp_GeVinv  = toml::find_or<double>(sTbl, "q_over_p_GeVinv", 1.0);
    double time_ns        = toml::find_or<double>(sTbl, "time_ns", 25.0);
    // Convert to Acts units:
    cfg.seeding.originStdDevLoc0   = loc0_mm   * 1_mm;
    cfg.seeding.originStdDevLoc1   = loc1_mm   * 1_mm;
    cfg.seeding.originStdDevPhi    = phi_deg   * 1_degree;
    cfg.seeding.originStdDevTheta  = theta_deg * 1_degree;
    cfg.seeding.originStdDevQOverP = qoverp_GeVinv / 1_GeV;
    cfg.seeding.originStdDevTime   = time_ns   * 1_ns;

    cfg.seeding.X0  = toml::find_or<double>(sTbl, "X0", 21.82);
    cfg.seeding.rho = toml::find_or<double>(sTbl, "rho", 2.329);
    cfg.seeding.x   = toml::find_or<double>(sTbl, "x", 25e-4);
    cfg.seeding.P   = toml::find_or<double>(sTbl, "P", 2500.);
    cfg.seeding.z   = toml::find_or<double>(sTbl, "z", 1.);

    cfg.seeding.nCellsThetaShort =
        toml::find_or<int>(sTbl, "n_cells_theta_short", 500);
    cfg.seeding.nCellsRhoShort =
        toml::find_or<int>(sTbl, "n_cells_rho_short", 4000);
    cfg.seeding.nCellsThetaLong =
        toml::find_or<int>(sTbl, "n_cells_theta_long", 500);
    cfg.seeding.nCellsRhoLong =
        toml::find_or<int>(sTbl, "n_cells_rho_long", 4000);
    cfg.seeding.nGLSIterations =
        toml::find_or<int>(sTbl, "n_gls_iterations", 2);
    cfg.seeding.minXCount =
        toml::find_or<int>(sTbl, "min_x_count", 3);
    cfg.seeding.minSeedSize =
        toml::find_or<int>(sTbl, "min_seed_size", 5);
    cfg.seeding.maxSeedSize =
        toml::find_or<int>(sTbl, "max_seed_size", 100);
    cfg.seeding.maxChi2 =
        toml::find_or<double>(sTbl, "max_chi2", 1e16);
  }

  // [TRACK_FITTING]
  if (root.contains("TRACK_FITTING")) {
    const auto& fTbl = toml::find(root, "TRACK_FITTING");
    cfg.fitting.maxSteps =
        toml::find_or<std::size_t>(fTbl, "max_steps", 100000u);
    cfg.fitting.resolvePassive =
        toml::find_or<bool>(fTbl, "resolve_passive", false);
    cfg.fitting.resolveMaterial =
        toml::find_or<bool>(fTbl, "resolve_material", true);
    cfg.fitting.resolveSensitive =
        toml::find_or<bool>(fTbl, "resolve_sensitive", true);
    cfg.fitting.inputCandidates =
        toml::find_or<std::string>(fTbl, "input_candidates", "Seeds");
    cfg.fitting.outputTracks =
        toml::find_or<std::string>(fTbl, "output_tracks", "Tracks");
  }

  // [MEASUREMENT_WRITER]
  if (root.contains("MEASUREMENT_WRITER")) {
    const auto& mTbl = toml::find(root, "MEASUREMENT_WRITER");
    cfg.measurementWriter.enable =
        toml::find_or<bool>(mTbl, "enable", false);
    cfg.measurementWriter.input =
        toml::find_or<std::string>(mTbl, "inputMeasurements", "Measurements");
    cfg.measurementWriter.treeName =
        toml::find_or<std::string>(mTbl, "treeName", "measurements");
    cfg.measurementWriter.outputFile =
        toml::find_or<std::string>(mTbl, "output_file", "measurements.root");
  }

  // [SEED_WRITER]
  if (root.contains("SEED_WRITER")) {
    const auto& sWTbl = toml::find(root, "SEED_WRITER");
    cfg.seedWriter.enable =
        toml::find_or<bool>(sWTbl, "enable", false);
    cfg.seedWriter.inputSeeds =
        toml::find_or<std::string>(sWTbl, "inputSeeds", "Seeds");
    cfg.seedWriter.inputTruthClusters =
        toml::find_or<std::string>(sWTbl, "inputTruthClusters", "SimClusters");
    cfg.seedWriter.treeName =
        toml::find_or<std::string>(sWTbl, "treeName", "seeds");
    cfg.seedWriter.outputFile =
        toml::find_or<std::string>(sWTbl, "output_file", "seeds-0.root");
  }

  // [TRACK_WRITER]
  if (root.contains("TRACK_WRITER")) {
    const auto& tWTbl = toml::find(root, "TRACK_WRITER");
    cfg.trackWriter.enable =
        toml::find_or<bool>(tWTbl, "enable", true);
    cfg.trackWriter.inputTracks =
        toml::find_or<std::string>(tWTbl, "inputTracks", "Tracks");
    cfg.trackWriter.inputSimClusters =
        toml::find_or<std::string>(tWTbl, "inputSimClusters", "SimClusters");
    cfg.trackWriter.treeName =
        toml::find_or<std::string>(tWTbl, "treeName", "fitted-tracks");
    cfg.trackWriter.outputFile =
        toml::find_or<std::string>(tWTbl, "output_file", "fitted-tracks-0.root");
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

int runFullSimTracking(const std::string& configPath) {
  FullSimTrackingRunConfig cfg;
  try {
    cfg = parseFullSimTrackingRunConfig(configPath);
  } catch (const std::exception& e) {
    std::cerr << "Error parsing FullSimTrackingRun config: " << e.what() << "\n";
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

  // Surface map for reader
  std::map<Acts::GeometryIdentifier, const Acts::Surface*> surfaceMap;
  for (const auto& vol : detector->volumes()) {
    for (const auto& surf : vol->surfaces()) {
      if (surf->geometryId().sensitive()) {
        surfaceMap[surf->geometryId()] = surf;
      }
    }
  }

  // Alignment store: from file or built
  std::shared_ptr<AlignmentContext::AlignmentStore> aStore;
  if (cfg.alignment.useFile && !cfg.alignment.filePath.empty()) {
    AlignmentParametersProvider::Config alignmentProviderCfg;
    alignmentProviderCfg.filePath = cfg.alignment.filePath;
    alignmentProviderCfg.treeName = cfg.alignment.treeName;
    AlignmentParametersProvider alignmentProvider(alignmentProviderCfg);
    aStore = alignmentProvider.getAlignmentStore();
  } else {
    aStore = detail::makeAlignmentStore(gctx, detector.get());
  }

  AlignmentContext alignCtx(aStore);
  gctx = Acts::GeometryContext{alignCtx};

  // Magnetic field
  auto field = E320Geometry::buildMagField(gctx);

  // Reference surfaces
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

  Acts::Transform3 seedingRefSurfTransform = Acts::Transform3::Identity();
  seedingRefSurfTransform.translation() =
      Acts::Vector3(goInst.ipTcDistance - 0.3_mm, 0, 0);
  seedingRefSurfTransform.rotate(refSurfToWorldRotationX);
  seedingRefSurfTransform.rotate(refSurfToWorldRotationY);
  seedingRefSurfTransform.rotate(refSurfToWorldRotationZ);

  auto seedingRefSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      seedingRefSurfTransform,
      std::make_shared<Acts::RectangleBounds>(halfX, halfY));

  Acts::GeometryIdentifier seedingRefSurfaceGeoId;
  seedingRefSurfaceGeoId.setExtra(1);
  seedingRefSurface->assignGeometryId(std::move(seedingRefSurfaceGeoId));

  Acts::Transform3 trackingRefSurfTransform = Acts::Transform3::Identity();
  trackingRefSurfTransform.translation() =
      Acts::Vector3(
          goInst.dipoleCenterPrimary + goInst.dipoleHalfPrimary + 0.01_mm, 0,
          0);
  trackingRefSurfTransform.rotate(refSurfToWorldRotationX);
  trackingRefSurfTransform.rotate(refSurfToWorldRotationY);
  trackingRefSurfTransform.rotate(refSurfToWorldRotationZ);

  auto trackingRefSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      trackingRefSurfTransform,
      std::make_shared<Acts::RectangleBounds>(halfX, halfY));

  Acts::GeometryIdentifier trackingRefSurfaceGeoId;
  trackingRefSurfaceGeoId.setExtra(2);
  trackingRefSurface->assignGeometryId(std::move(trackingRefSurfaceGeoId));

  // --------------------------------------------------------------
  // Event reading
  SimpleSourceLink::SurfaceAccessor surfaceAccessor{detector.get()};

  Sequencer::Config seqCfg;
  seqCfg.numThreads = cfg.sequencer.numThreads;
  seqCfg.skip       = cfg.sequencer.skip;
  seqCfg.trackFpes  = cfg.sequencer.trackFpes;
  seqCfg.logLevel   = logLevel;
  Sequencer sequencer(seqCfg);

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));

  // Sim data reader
  RootSimClusterReader::Config readerCfg;
  readerCfg.outputSourceLinks    = cfg.reader.outputSourceLinks;
  readerCfg.outputSimClusters    = cfg.reader.outputSimClusters;
  readerCfg.treeName             = cfg.reader.treeName;
  readerCfg.minGeoId             = cfg.reader.minGeoId;
  readerCfg.maxGeoId             = cfg.reader.maxGeoId;
  readerCfg.surfaceLocalToGlobal = cfg.reader.surfaceLocalToGlobal;
  readerCfg.surfaceMap           = surfaceMap;

  for (const auto& entry :
       std::filesystem::directory_iterator(cfg.reader.clustersDir)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".root") {
      continue;
    }
    readerCfg.filePaths.push_back(entry.path().string());
  }

  sequencer.addReader(
      std::make_shared<RootSimClusterReader>(readerCfg, logLevel));

  // --------------------------------------------------------------
  // HT seeding setup
  double X0  = cfg.seeding.X0;
  double rho = cfg.seeding.rho;
  double x   = cfg.seeding.x;
  double P   = cfg.seeding.P;
  double z   = cfg.seeding.z;

  double t        = rho * x;  // [g / cm2]
  double thetaRms = 13.6 / P * z * std::sqrt(t / X0) *
                    (1 + 0.038 * std::log(t * z * z / (X0)));

  HoughTransformSeeder::Config htSeederCfg;
  htSeederCfg.boundBoxHalfPrimary = goInst.tcHalfPrimary;
  htSeederCfg.boundBoxHalfLong    = goInst.tcHalfLong;
  htSeederCfg.boundBoxHalfShort   = goInst.tcHalfShort;

  htSeederCfg.nCellsThetaShort = cfg.seeding.nCellsThetaShort;
  htSeederCfg.nCellsRhoShort   = cfg.seeding.nCellsRhoShort;
  htSeederCfg.nCellsThetaLong  = cfg.seeding.nCellsThetaLong;
  htSeederCfg.nCellsRhoLong    = cfg.seeding.nCellsRhoLong;
  htSeederCfg.nGLSIterations   = cfg.seeding.nGLSIterations;

  htSeederCfg.primaryIdx = goInst.primaryIdx;
  htSeederCfg.longIdx    = goInst.longIdx;
  htSeederCfg.shortIdx   = goInst.shortIdx;

  HoughTransformSeeder::Options htSeederOpt;
  htSeederOpt.boundBoxCenterPrimary = goInst.tcCenterPrimary;
  htSeederOpt.boundBoxCenterLong    = goInst.tcCenterLong;
  htSeederOpt.boundBoxCenterShort   = goInst.tcCenterShort;

  htSeederOpt.firstLayerId = goInst.tcParameters.front().geoId;
  htSeederOpt.lastLayerId  = goInst.tcParameters.back().geoId;
  htSeederOpt.nLayers      = goInst.tcParameters.size();

  htSeederOpt.minXCount   = cfg.seeding.minXCount;
  htSeederOpt.minSeedSize = cfg.seeding.minSeedSize;
  htSeederOpt.maxSeedSize = cfg.seeding.maxSeedSize;

  htSeederOpt.primaryInterchipDistance = goInst.interChipDistance;
  htSeederOpt.thetaRms                 = thetaRms;
  htSeederOpt.maxChi2                  = cfg.seeding.maxChi2;

  E320SeedingAlgorithm::Config seedingAlgoCfg;
  seedingAlgoCfg.htSeeder = std::make_shared<HoughTransformSeeder>(htSeederCfg);
  seedingAlgoCfg.htOptions        = htSeederOpt;
  seedingAlgoCfg.inputSourceLinks = cfg.seeding.inputSourceLinks;
  seedingAlgoCfg.outputSeeds      = cfg.seeding.outputSeeds;
  seedingAlgoCfg.minLayers        = cfg.seeding.minLayers;
  seedingAlgoCfg.maxLayers        = cfg.seeding.maxLayers;
  seedingAlgoCfg.beamlineTilt     = cfg.seeding.beamlineTilt;
  seedingAlgoCfg.referenceSurface = seedingRefSurface.get();

  Acts::BoundVector trackOriginStdDevPrior;
  trackOriginStdDevPrior[Acts::eBoundLoc0]   = cfg.seeding.originStdDevLoc0;
  trackOriginStdDevPrior[Acts::eBoundLoc1]   = cfg.seeding.originStdDevLoc1;
  trackOriginStdDevPrior[Acts::eBoundTime]   = cfg.seeding.originStdDevTime;
  trackOriginStdDevPrior[Acts::eBoundPhi]    = cfg.seeding.originStdDevPhi;
  trackOriginStdDevPrior[Acts::eBoundTheta]  = cfg.seeding.originStdDevTheta;
  trackOriginStdDevPrior[Acts::eBoundQOverP] = cfg.seeding.originStdDevQOverP;
  seedingAlgoCfg.originCov =
      trackOriginStdDevPrior.cwiseProduct(trackOriginStdDevPrior).asDiagonal();

  sequencer.addAlgorithm(
      std::make_shared<E320SeedingAlgorithm>(seedingAlgoCfg, logLevel));

  // --------------------------------------------------------------
  // Track fitting
  Acts::GainMatrixUpdater  kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;

  Acts::KalmanFitterExtensions<RecoTrajectory> extensions;
  extensions.calibrator.connect<&simpleSourceLinkCalibrator<RecoTrajectory>>();
  extensions.updater
      .connect<&Acts::GainMatrixUpdater::operator()<RecoTrajectory>>(
          &kfUpdater);
  extensions.smoother
      .connect<&Acts::GainMatrixSmoother::operator()<RecoTrajectory>>(
          &kfSmoother);
  extensions.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  auto propOptions = PropagatorOptions(gctx, mctx);
  propOptions.maxSteps = cfg.fitting.maxSteps;

  auto options =
      Acts::KalmanFitterOptions(gctx, mctx, cctx, extensions, propOptions);
  options.referenceSurface = trackingRefSurface.get();

  Acts::Experimental::DetectorNavigator::Config navCfg;
  navCfg.detector         = detector.get();
  navCfg.resolvePassive   = cfg.fitting.resolvePassive;
  navCfg.resolveMaterial  = cfg.fitting.resolveMaterial;
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
      .kfOptions            = options};

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

  // Optional measurement writer (usually off for sim)
  if (cfg.measurementWriter.enable) {
    auto measurementWriterCfg = RootMeasurementWriter::Config();
    measurementWriterCfg.inputMeasurements = cfg.measurementWriter.input;
    measurementWriterCfg.treeName          = cfg.measurementWriter.treeName;
    measurementWriterCfg.filePath =
        (outDir / cfg.measurementWriter.outputFile).string();
    sequencer.addWriter(
        std::make_shared<RootMeasurementWriter>(measurementWriterCfg, logLevel));
  }

  // Seed writer
  if (cfg.seedWriter.enable) {
    auto seedWriterCfg = RootSimSeedWriter::Config();
    seedWriterCfg.inputSeeds         = cfg.seedWriter.inputSeeds;
    seedWriterCfg.inputTruthClusters = cfg.seedWriter.inputTruthClusters;
    seedWriterCfg.treeName           = cfg.seedWriter.treeName;
    seedWriterCfg.filePath           =
        (outDir / cfg.seedWriter.outputFile).string();

    sequencer.addWriter(
        std::make_shared<RootSimSeedWriter>(seedWriterCfg, logLevel));
  }

  // Track writer
  if (cfg.trackWriter.enable) {
    auto trackWriterCfg = RootSimTrackWriter::Config();
    trackWriterCfg.surfaceAccessor
        .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
            &surfaceAccessor);
    trackWriterCfg.referenceSurface = trackingRefSurface.get();
    trackWriterCfg.inputTracks      = cfg.trackWriter.inputTracks;
    trackWriterCfg.inputSimClusters = cfg.trackWriter.inputSimClusters;
    trackWriterCfg.treeName         = cfg.trackWriter.treeName;
    trackWriterCfg.filePath         =
        (outDir / cfg.trackWriter.outputFile).string();

    sequencer.addWriter(
        std::make_shared<RootSimTrackWriter>(trackWriterCfg, logLevel));
  }

  return sequencer.run();
}

} // namespace TrackingPipeline

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <FullSimTrackingRun.conf>\n";
    return 1;
  }
  std::string confPath = argv[1];
  return TrackingPipeline::runFullSimTracking(confPath);
}
