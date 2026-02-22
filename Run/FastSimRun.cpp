#include "TrackingPipeline/Run/FastSimRun.hpp"

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/RectangleBounds.hpp>

#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <toml.hpp>

#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Alignment/detail/AlignmentStoreBuilders.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/AlignmentParametersProvider.hpp"
#include "TrackingPipeline/Io/DummyReader.hpp"
#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"
#include "TrackingPipeline/Simulation/GaussianVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"
#include "TrackingPipeline/Simulation/SimpleDigitizer.hpp"
#include "TrackingPipeline/Simulation/SphericalMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/UniformBackgroundCreator.hpp"

using namespace Acts::UnitLiterals;
namespace eg = E320Geometry;

// Propagator short-hands
using ActionList = Acts::ActionList<>;
using AbortList  = Acts::AbortList<Acts::EndOfWorldReached>;

using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                    Acts::Experimental::DetectorNavigator>;
using PropagatorOptions =
    typename Propagator::template Options<ActionList, AbortList>;

using RecoTrajectory     = Acts::VectorMultiTrajectory;
using RecoTrackContainer = Acts::VectorTrackContainer;
using KF                 = Acts::KalmanFitter<Propagator, RecoTrajectory>;

std::unique_ptr<const eg::GeometryOptions> eg::GeometryOptions::m_instance =
    nullptr;

namespace TrackingPipeline {

FastSimRunConfig parseFastSimRunConfig(const std::string& path) {
  auto root = toml::parse(path);
  FastSimRunConfig cfg;
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

  // [DUMMY_READER]
  if (root.contains("DUMMY_READER")) {
    const auto& dTbl = toml::find(root, "DUMMY_READER");
    cfg.dummyReader.nEvents =
        toml::find_or<std::size_t>(dTbl, "n_events", 100000u);
    cfg.dummyReader.outputSourceLinks =
        toml::find_or<std::string>(dTbl, "output_source_links", "SimMeasurements");
    cfg.dummyReader.outputSimClusters =
        toml::find_or<std::string>(dTbl, "output_sim_clusters", "SimClusters");
  }

  // [SIGNAL]
  if (root.contains("SIGNAL")) {
    const auto& sTbl = toml::find(root, "SIGNAL");
    cfg.signal.resolvePassive =
        toml::find_or<bool>(sTbl, "resolve_passive", false);
    cfg.signal.resolveMaterial =
        toml::find_or<bool>(sTbl, "resolve_material", true);
    cfg.signal.resolveSensitive =
        toml::find_or<bool>(sTbl, "resolve_sensitive", true);
    cfg.signal.resY_um =
        toml::find_or<double>(sTbl, "res_y_um", 5.0);
    cfg.signal.resZ_um =
        toml::find_or<double>(sTbl, "res_z_um", 5.0);
    cfg.signal.pMin_GeV =
        toml::find_or<double>(sTbl, "p_min_GeV", 2.0);
    cfg.signal.pMax_GeV =
        toml::find_or<double>(sTbl, "p_max_GeV", 3.0);
    cfg.signal.phiMin =
        toml::find_or<double>(sTbl, "phi_min", -0.001);
    cfg.signal.phiMax =
        toml::find_or<double>(sTbl, "phi_max",  0.001);
    cfg.signal.thetaMin =
        toml::find_or<double>(sTbl, "theta_min", 1.57079632679 - 0.003);
    cfg.signal.thetaMax =
        toml::find_or<double>(sTbl, "theta_max", 1.57079632679 + 0.003);
    cfg.signal.vtxSigma_um =
        toml::find_or<double>(sTbl, "vtx_sigma_um", 30.0);
  }

  // [SIGNAL_MEAS_CREATOR]
  if (root.contains("SIGNAL_MEAS_CREATOR")) {
    const auto& mTbl = toml::find(root, "SIGNAL_MEAS_CREATOR");
    cfg.signalMeasCreator.inputSourceLinks =
        toml::find_or<std::string>(mTbl, "input_source_links", "SimMeasurements");
    cfg.signalMeasCreator.inputSimClusters =
        toml::find_or<std::string>(mTbl, "input_sim_clusters", "SimClusters");
    cfg.signalMeasCreator.outputSourceLinks =
        toml::find_or<std::string>(mTbl, "output_source_links", "Measurements");
    cfg.signalMeasCreator.outputSimClusters =
        toml::find_or<std::string>(mTbl, "output_sim_clusters", "Clusters");
    cfg.signalMeasCreator.maxSteps =
        toml::find_or<std::size_t>(mTbl, "max_steps", 1000u);
    cfg.signalMeasCreator.nMeasurements =
        toml::find_or<int>(mTbl, "n_measurements", 1);
    cfg.signalMeasCreator.minSensitiveId =
        toml::find_or<int>(mTbl, "min_sensitive_id", 40);
  }

  // [BACKGROUND]
  if (root.contains("BACKGROUND")) {
    const auto& bTbl = toml::find(root, "BACKGROUND");
    cfg.background.enable =
        toml::find_or<bool>(bTbl, "enable", false);
    cfg.background.nHits =
        toml::find_or<std::size_t>(bTbl, "n_hits", 700u);
    cfg.background.resY_um =
        toml::find_or<double>(bTbl, "res_y_um", 5.0);
    cfg.background.resZ_um =
        toml::find_or<double>(bTbl, "res_z_um", 5.0);
    cfg.background.inputSourceLinks =
        toml::find_or<std::string>(bTbl, "input_source_links", "Measurements1");
    cfg.background.inputSimClusters =
        toml::find_or<std::string>(bTbl, "input_sim_clusters", "Clusters1");
    cfg.background.outputSourceLinks =
        toml::find_or<std::string>(bTbl, "output_source_links", "Measurements");
    cfg.background.outputSimClusters =
        toml::find_or<std::string>(bTbl, "output_sim_clusters", "Clusters");
    cfg.background.nMeasurements =
        toml::find_or<int>(bTbl, "n_measurements", 1);
  }

  // [CLUSTER_WRITER]
  if (root.contains("CLUSTER_WRITER")) {
    const auto& wTbl = toml::find(root, "CLUSTER_WRITER");
    cfg.clusterWriter.enable =
        toml::find_or<bool>(wTbl, "enable", true);
    cfg.clusterWriter.inputClusters =
        toml::find_or<std::string>(wTbl, "inputClusters", "Clusters");
    cfg.clusterWriter.treeName =
        toml::find_or<std::string>(wTbl, "treeName", "clusters");
    cfg.clusterWriter.filePrefix =
        toml::find_or<std::string>(wTbl, "file_prefix", "clusters-");
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

int runFastSim(const std::string& configPath, int jobId) {
  FastSimRunConfig cfg;
  try {
    cfg = parseFastSimRunConfig(configPath);
  } catch (const std::exception& e) {
    std::cerr << "Error parsing FastSimRun config: " << e.what() << "\n";
    return 1;
  }

  Acts::Logging::Level logLevel = getLogLevel(cfg.main.logLevel);
  const auto& goInst            = *eg::GeometryOptions::instance();

  // Contexts
  Acts::GeometryContext      gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext   cctx;

  // Detector
  auto detector = eg::buildDetector(gctx);

  std::vector<const Acts::Surface*> detSurfaces;
  for (const auto* vol : detector->volumes()) {
    for (const auto* surf : vol->surfaces()) {
      if (surf->geometryId().sensitive()) {
        detSurfaces.push_back(surf);
      }
    }
  }

  // Alignment: either from file or Gaussian pattern as before
  std::shared_ptr<AlignmentContext::AlignmentStore> aStore;
  if (cfg.alignment.useFile && !cfg.alignment.filePath.empty()) {
    AlignmentParametersProvider::Config alignmentProviderCfg;
    alignmentProviderCfg.filePath = cfg.alignment.filePath;
    alignmentProviderCfg.treeName = cfg.alignment.treeName;
    AlignmentParametersProvider alignmentProvider(alignmentProviderCfg);
    aStore = alignmentProvider.getAlignmentStore();
  } else {
    Acts::Vector3 globalShiftMean(0, 0_mm, 0_mm);
    Acts::Vector3 globalShiftStdErr(0, 0_mm, 0_mm);

    std::unordered_map<int, Acts::Vector3> localShiftsMean{
        {10, Acts::Vector3(0_mm, 0_um, 0_um)},
        {12, Acts::Vector3(0_mm, 0_um, 0_um)},
        {14, Acts::Vector3(0_mm, 0_um, 0_um)},
        {16, Acts::Vector3(0_mm, 0_um, 0_um)},
        {18, Acts::Vector3(0_mm, 0_um, 0_um)}};
    std::unordered_map<int, Acts::Vector3> localShiftsStdErr{
        {10, Acts::Vector3(0_mm, 30_um, 30_um)},
        {12, Acts::Vector3(0_mm, 30_um, 30_um)},
        {14, Acts::Vector3(0_mm, 30_um, 30_um)},
        {16, Acts::Vector3(0_mm, 30_um, 30_um)},
        {18, Acts::Vector3(0_mm, 30_um, 30_um)}};

    Acts::Vector3 globalAnglesMean(0_rad, 0_rad, 0_rad);
    Acts::Vector3 globalAnglesStdErr(0_rad, 0_rad, 0_rad);

    std::unordered_map<int, Acts::Vector3> localAnglesMean{
        {10, Acts::Vector3(0_rad, 0_rad, 0_rad)},
        {12, Acts::Vector3(0_rad, 0_rad, 0_rad)},
        {14, Acts::Vector3(0_rad, 0_rad, 0_rad)},
        {16, Acts::Vector3(0_rad, 0_rad, 0_rad)},
        {18, Acts::Vector3(0_rad, 0_rad, 0_rad)}};
    std::unordered_map<int, Acts::Vector3> localAnglesStdErr{
        {10, Acts::Vector3(0_rad, 0_rad, 1e-3_rad)},
        {12, Acts::Vector3(0_rad, 0_rad, 1e-3_rad)},
        {14, Acts::Vector3(0_rad, 0_rad, 1e-3_rad)},
        {16, Acts::Vector3(0_rad, 0_rad, 1e-3_rad)},
        {18, Acts::Vector3(0_rad, 0_rad, 1e-3_rad)}};

    aStore = detail::makeAlignmentStore(
        gctx, detector.get(), globalShiftMean, globalAnglesStdErr,
        localShiftsMean, localShiftsStdErr, globalAnglesMean, globalAnglesStdErr,
        localAnglesMean, localAnglesStdErr);
  }

  AlignmentContext alignCtx(aStore);
  gctx = Acts::GeometryContext{alignCtx};

  // Magnetic field
  auto field = eg::buildMagField(gctx);

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

  Acts::Transform3 refSurfTransform = Acts::Transform3::Identity();
  refSurfTransform.translation() = Acts::Vector3::Zero();
  refSurfTransform.rotate(refSurfToWorldRotationX);
  refSurfTransform.rotate(refSurfToWorldRotationY);
  refSurfTransform.rotate(refSurfToWorldRotationZ);

  auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      refSurfTransform, std::make_shared<Acts::RectangleBounds>(halfX, halfY));

  Acts::GeometryIdentifier geoId;
  geoId.setExtra(1);
  refSurface->assignGeometryId(std::move(geoId));

  // --------------------------------------------------------------
  // Event reading
  SimpleSourceLink::SurfaceAccessor surfaceAccessor{detector.get()};

  // Sequencer
  Sequencer::Config seqCfg;
  seqCfg.numThreads = cfg.sequencer.numThreads;
  seqCfg.skip       = cfg.sequencer.skip + jobId * cfg.dummyReader.nEvents;
  seqCfg.trackFpes  = cfg.sequencer.trackFpes;
  seqCfg.logLevel   = logLevel;
  Sequencer sequencer(seqCfg);

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));

  // Dummy reader
  DummyReader::Config dummyReaderCfg;
  dummyReaderCfg.outputSourceLinks = cfg.dummyReader.outputSourceLinks;
  dummyReaderCfg.outputSimClusters = cfg.dummyReader.outputSimClusters;
  dummyReaderCfg.nEvents           = (jobId + 1) * cfg.dummyReader.nEvents;

  sequencer.addReader(std::make_shared<DummyReader>(dummyReaderCfg));

  // --------------------------------------------------------------
  // Measurements creator (signal)
  Acts::Experimental::DetectorNavigator::Config cptNavCfg;
  cptNavCfg.detector        = detector.get();
  cptNavCfg.resolvePassive  = cfg.signal.resolvePassive;
  cptNavCfg.resolveMaterial = cfg.signal.resolveMaterial;
  cptNavCfg.resolveSensitive = cfg.signal.resolveSensitive;

  Acts::Experimental::DetectorNavigator measCreatorNavigator(
      cptNavCfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));
  Acts::EigenStepper<> measCreatorStepper(field);

  Propagator measCreatorPropagator(std::move(measCreatorStepper),
                                   std::move(measCreatorNavigator));

  // Digitizer
  SimpleDigitizer::Config digitizerCfg;
  digitizerCfg.resolution = {cfg.signal.resY_um * 1_um, cfg.signal.resZ_um * 1_um};
  auto digitizer = std::make_shared<SimpleDigitizer>(digitizerCfg);

  // Vertex generator
  GaussianVertexGenerator::Config vertexGenCfg;
  vertexGenCfg.mean = Acts::Vector3(0, 0, 0);
  vertexGenCfg.cov  = Acts::SquareMatrix3::Identity() *
                      (cfg.signal.vtxSigma_um * 1_um);
  auto vertexGen = std::make_shared<GaussianVertexGenerator>(vertexGenCfg);

  // Momentum generator
  SphericalMomentumGenerator::Config momGenCfg;
  momGenCfg.pRange     = {cfg.signal.pMin_GeV * 1_GeV,
                          cfg.signal.pMax_GeV * 1_GeV};
  momGenCfg.phiRange   = {cfg.signal.phiMin, cfg.signal.phiMax};
  momGenCfg.thetaRange = {cfg.signal.thetaMin, cfg.signal.thetaMax};

  auto momGen = std::make_shared<SphericalMomentumGenerator>(momGenCfg);

  // ----------------------------------------------
  // Measurement creator (signal)
  MeasurementsCreator::Config measCreatorCfg{
      .vertexGenerator   = vertexGen,
      .momentumGenerator = momGen,
      .hitDigitizer      = digitizer,
      .referenceSurface  = refSurface.get(),
      .maxSteps          = static_cast<std::size_t>(cfg.signalMeasCreator.maxSteps),
      .isSignal          = true,
      .hypothesis        = Acts::ParticleHypothesis::electron(),
      .charge            = -1_e};

  std::unordered_map<Acts::GeometryIdentifier, MeasurementsCreator::Constraints>
      measCreatorConstraints;
  for (auto& det : detector->detectorElements()) {
    const auto& surface = det->surface();
    const auto& gid     = surface.geometryId();
    if (gid.sensitive() && gid.sensitive() >= cfg.signalMeasCreator.minSensitiveId) {
      measCreatorConstraints.insert({gid, {-3, 3, -3, 3}});
    }
  }
  measCreatorCfg.constraints = measCreatorConstraints;

  auto measCreator = std::make_shared<MeasurementsCreator>(
      measCreatorPropagator, measCreatorCfg);

  MeasurementsEmbeddingAlgorithm::Config measCreatorAlgoCfg;
  measCreatorAlgoCfg.inputSourceLinks  = cfg.signalMeasCreator.inputSourceLinks;
  measCreatorAlgoCfg.inputSimClusters  = cfg.signalMeasCreator.inputSimClusters;
  measCreatorAlgoCfg.outputSourceLinks = cfg.signalMeasCreator.outputSourceLinks;
  measCreatorAlgoCfg.outputSimClusters = cfg.signalMeasCreator.outputSimClusters;
  measCreatorAlgoCfg.measurementGenerator = measCreator;
  measCreatorAlgoCfg.randomNumberSvc =
      std::make_shared<RandomNumbers>(RandomNumbers::Config());
  measCreatorAlgoCfg.nMeasurements = cfg.signalMeasCreator.nMeasurements;

  sequencer.addAlgorithm(std::make_shared<MeasurementsEmbeddingAlgorithm>(
      measCreatorAlgoCfg, logLevel));

  // --------------------------------------------------------------
  // Background (optional)
  if (cfg.background.enable) {
    UniformBackgroundCreator::Config bkgCreatorCfg;
    bkgCreatorCfg.resolution = {cfg.background.resY_um * 1_um,
                                cfg.background.resZ_um * 1_um};
    bkgCreatorCfg.nMeasurements = cfg.background.nHits;
    bkgCreatorCfg.surfaces      = detSurfaces;

    auto bkgCreator = std::make_shared<UniformBackgroundCreator>(bkgCreatorCfg);

    MeasurementsEmbeddingAlgorithm::Config bkgCreatorAlgoCfg;
    bkgCreatorAlgoCfg.inputSourceLinks  = cfg.background.inputSourceLinks;
    bkgCreatorAlgoCfg.inputSimClusters  = cfg.background.inputSimClusters;
    bkgCreatorAlgoCfg.outputSourceLinks = cfg.background.outputSourceLinks;
    bkgCreatorAlgoCfg.outputSimClusters = cfg.background.outputSimClusters;
    bkgCreatorAlgoCfg.measurementGenerator = bkgCreator;
    bkgCreatorAlgoCfg.randomNumberSvc =
        std::make_shared<RandomNumbers>(RandomNumbers::Config());
    bkgCreatorAlgoCfg.nMeasurements = cfg.background.nMeasurements;

    // Add if/when you want background enabled
    // sequencer.addAlgorithm(std::make_shared<MeasurementsEmbeddingAlgorithm>(
    //     bkgCreatorAlgoCfg, logLevel));
  }

  // Output dir
  std::filesystem::path outDir(cfg.main.outputDir);
  try {
    std::filesystem::create_directories(outDir);
  } catch (const std::exception& e) {
    std::cerr << "Error creating output directory '" << outDir.string()
              << "': " << e.what() << "\n";
    return 1;
  }

  // --------------------------------------------------------------
  // Sim cluster writer
  if (cfg.clusterWriter.enable) {
    auto clusterWriterCfgSig = RootSimClusterWriter::Config();
    clusterWriterCfgSig.inputClusters = cfg.clusterWriter.inputClusters;
    clusterWriterCfgSig.treeName      = cfg.clusterWriter.treeName;
    clusterWriterCfgSig.filePath      =
        (outDir / (cfg.clusterWriter.filePrefix +
                    std::to_string(jobId) + ".root")).string();

    sequencer.addWriter(
        std::make_shared<RootSimClusterWriter>(clusterWriterCfgSig, logLevel));
    }

  return sequencer.run();
}

} // namespace TrackingPipeline

int main(int argc, char* argv[]) {
  std::string confPath;
  int jobId = 0;

  if (argc == 2) {
    // Only config path given
    confPath = argv[1];
  } else if (argc >= 3) {
    // jobId + config
    jobId    = std::stoi(argv[1]);
    confPath = argv[2];
  } else {
    std::cerr << "Usage: " << argv[0]
              << " <FastSimRun.conf> [job_id]\n";
    return 1;
  }

  return TrackingPipeline::runFastSim(confPath, jobId);
}

