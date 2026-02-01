#include "TrackingPipeline/Run/PipelineRun.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
// Geometry & Infrastructure
#include "TrackingPipeline/Geometry/ApollonGeometry.hpp"
#include "TrackingPipeline/Geometry/ApollonGeometryConstraints.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/ReaderRegistry.hpp"
#include "TrackingPipeline/Infrastructure/WriterRegistry.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmRegistry.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Alignment/detail/AlignmentUtils.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Io/AlignmentParametersProvider.hpp"
#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
//preprocessing
#include "TrackingPipeline/Preprocessing/Preprocessing.hpp"
// Fitting services (singleton with propagator + KF)
#include "TrackingPipeline/TrackFitting/FittingServices.hpp"
#include "TrackingPipeline/TrackFitting/KFTrackFittingAlgorithm.hpp"

// Standard includes
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <toml.hpp>

namespace TrackingPipeline {

using namespace Acts::UnitLiterals;
using detail::makeAlignmentStore;

//-------------------- helpers --------------------
static Acts::BinningValue parseBinningValue(const std::string& name,
                                            Acts::BinningValue fallback) {
  if (name == "binX") return Acts::BinningValue::binX;
  if (name == "binY") return Acts::BinningValue::binY;
  if (name == "binZ") return Acts::BinningValue::binZ;
  // extend here if you ever need more
  return fallback;
}

//-------------------- parsing config file --------------------
PipelineConfig parsePipelineConfig(const std::string& path) {
  auto root = toml::parse(path);
  PipelineConfig cfg;
	cfg.root = root;
	//[MAIN]
  cfg.main.outputDirLoc  = toml::find<std::string>(root, "MAIN", "output_dir_loc");
  cfg.main.outputDirName = toml::find<std::string>(root, "MAIN", "output_dir_name");
  cfg.main.logLevel      = toml::find<std::string>(root, "MAIN", "log_level");
	//[DETECTOR]
	const auto& detTbl = toml::find(root, "DETECTOR");
  cfg.detector.fillSurfaceMap =
      toml::find<bool>(detTbl, "fill_surface_map");
	cfg.detector.storeMode =
      toml::find_or<std::string>(detTbl, "store_mode", std::string("fixed"));
  cfg.detector.longBinningName =
      toml::find_or<std::string>(detTbl, "long_binning", std::string{});
  cfg.detector.shortBinningName =
      toml::find_or<std::string>(detTbl, "short_binning", std::string{});
  cfg.detector.longTransStdUm =
      toml::find_or<double>(detTbl, "long_trans_std", 0.0);
  cfg.detector.shortTransStdUm =
      toml::find_or<double>(detTbl, "short_trans_std", 0.0);

	// [ALIGNMENT_PROVIDER] (optional)
  if (root.contains("ALIGNMENT_PROVIDER")) {
    const auto& provTbl = toml::find(root, "ALIGNMENT_PROVIDER");
    cfg.alignmentProvider.enable =
        toml::find_or<bool>(provTbl, "enable", false);
    cfg.alignmentProvider.treeName =
        toml::find_or<std::string>(provTbl, "tree_name", std::string("alignment-parameters"));
    cfg.alignmentProvider.parametersFilepath =
        toml::find_or<std::string>(provTbl, "parameters_filepath", std::string{});
  } else {
    cfg.alignmentProvider.enable = false;
    cfg.alignmentProvider.treeName.clear();
    cfg.alignmentProvider.parametersFilepath.clear();
  }

	//[DATA_READER]
  cfg.dataReader.type =
      toml::find<std::string>(root, "DATA_READER", "type");

	//[DATA_WRITER]
  const auto& writerTbl = toml::find(root, "DATA_WRITER");
  cfg.dataWriter.types  =
      toml::find<std::vector<std::string>>(writerTbl, "types");

  //[SEEDING]
  if (root.contains("SEEDING")) {
    const auto& seedingTbl = toml::find(root, "SEEDING");
    cfg.seeding.enable =
        toml::find_or<bool>(seedingTbl, "enable", false);
    cfg.seeding.types =
        toml::find_or<std::vector<std::string>>(seedingTbl, "types", {});
  } else {
    cfg.seeding.enable = false;
    cfg.seeding.types.clear();
  }

  //[TRACK_FITTING]]
  if (root.contains("TRACK_FITTING")) {
    const auto& fitTbl = toml::find(root, "TRACK_FITTING");
    cfg.trackFitting.enable =
        toml::find_or<bool>(fitTbl, "enable", false);
    cfg.trackFitting.types =
        toml::find_or<std::vector<std::string>>(fitTbl, "types", {});
  } else {
    cfg.trackFitting.enable = false;
    cfg.trackFitting.types.clear();
  }

  // [TRACK_CLEANING]
  if (root.contains("TRACK_CLEANING")) {
    const auto& clnTbl = toml::find(root, "TRACK_CLEANING");
    cfg.trackCleaning.enable =
        toml::find_or<bool>(clnTbl, "enable", false);
    cfg.trackCleaning.types =
        toml::find_or<std::vector<std::string>>(clnTbl, "types", {});
  } else {
    cfg.trackCleaning.enable = false;
    cfg.trackCleaning.types.clear();
  }

  // [PREPROCESSING] (optional)
  if (root.contains("PREPROCESSING")) {
    const auto& preTbl = toml::find(root, "PREPROCESSING");
    cfg.preprocessing.enable =
        toml::find_or<bool>(preTbl, "enable", false);
    cfg.preprocessing.inputDirs =
        toml::find_or<std::vector<std::string>>(preTbl, "input_dirs", {});
    cfg.preprocessing.inputTreeName =
        toml::find_or<std::string>(preTbl, "input_tree", std::string{});
    cfg.preprocessing.inputBranchName =
        toml::find_or<std::string>(preTbl, "input_branch", std::string("event"));
    cfg.preprocessing.outputFile =
        toml::find_or<std::string>(preTbl, "output_file", std::string("preprocessed.root"));
    cfg.preprocessing.outputTreeName =
        toml::find_or<std::string>(preTbl, "output_tree", std::string("MyTree"));
    cfg.preprocessing.skipEntries =
        toml::find_or<std::size_t>(preTbl, "skip_entries", 0u);
  } else {
    cfg.preprocessing.enable = false;
  }

  // [MEAS_EMBEDDING] (optional)
  if (root.contains("MEAS_EMBEDDING")) {
    const auto& measTbl = toml::find(root, "MEAS_EMBEDDING");
    cfg.measEmbedding.enable =
        toml::find_or(measTbl, "enable", false);
    cfg.measEmbedding.types =
        toml::find_or<std::vector<std::string>>(measTbl, "types", {});
  } else {
    cfg.measEmbedding.enable = false;
    cfg.measEmbedding.types.clear();
  }

  return cfg;
}
// Map string log level to Acts::Logging::Level
Acts::Logging::Level getLogLevel(const std::string& levelStr) {
  if (levelStr == "VERBOSE") return Acts::Logging::VERBOSE;
  if (levelStr == "DEBUG")   return Acts::Logging::DEBUG;
  if (levelStr == "INFO")    return Acts::Logging::INFO;
  if (levelStr == "WARNING") return Acts::Logging::WARNING;
  if (levelStr == "ERROR")   return Acts::Logging::ERROR;
  if (levelStr == "FATAL")   return Acts::Logging::FATAL;
  return Acts::Logging::INFO;
}

//-------------------- run pipeline --------------------
int runPipeline(const std::string& configPath) {
  // Parse configuration
  PipelineConfig pcfg;
  try {
    pcfg = parsePipelineConfig(configPath);
  } catch (const std::exception& e) {
    std::cerr << "Error parsing config: " << e.what() << "\n";
    return 1;
  }
	// Set global log level
  Acts::Logging::Level globalLogLevel = getLogLevel(pcfg.main.logLevel);

  // Compute run root from [MAIN] and ensure it exists
  std::string runRoot = pcfg.main.outputDirLoc + "/" + pcfg.main.outputDirName;
  try {
    std::filesystem::create_directories(runRoot);
  } catch (const std::exception& e) {
    std::cerr << "Error creating run directory '" << runRoot
              << "': " << e.what() << "\n";
    return 1;
  }

  // Contexts
  Acts::GeometryContext      gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext   cctx;

  // Detector setup
  auto detector = ApollonGeometry::buildDetector(gctx);

  SurfaceMap surfaceMap;
  if (pcfg.detector.fillSurfaceMap) {
    for (const auto& vol : detector->volumes()) {
      for (const auto& surf : vol->surfaces()) {
        if (surf->geometryId().sensitive()) {
          surfaceMap[surf->geometryId()] = surf;
        }
      }
    }
  }

	//--------------- Alignment store setup ---------------
	std::shared_ptr<AlignmentContext::AlignmentStore> aStore;
  const auto& goInst = *ApollonGeometry::GeometryOptions::instance();

  if (pcfg.detector.storeMode == "gaussian") {
    Acts::BinningValue longBinValue =
        parseBinningValue(pcfg.detector.longBinningName,  goInst.longBinValue);
    Acts::BinningValue shortBinValue =
        parseBinningValue(pcfg.detector.shortBinningName, goInst.shortBinValue);

    std::size_t longIdx  = detail::binningValueToIndex(longBinValue);
    std::size_t shortIdx = detail::binningValueToIndex(shortBinValue);
    // cfg values are in um, convert to acts mm
    double longTransStd  = pcfg.detector.longTransStdUm  * 1_um;
    double shortTransStd = pcfg.detector.shortTransStdUm * 1_um;
    aStore = 
			makeAlignmentStore(detector.get(), longIdx, longTransStd, shortIdx, shortTransStd);

  } else if (pcfg.detector.storeMode == "none") {
    // no internal misalignment pattern at all
    aStore = std::make_shared<AlignmentContext::AlignmentStore>();

  } else {
    aStore = makeAlignmentStore(detector.get()); // "fixed" pattern as in FastSimRun
  }

	//--------------- Alignment parameters provider ---------------
	if (pcfg.alignmentProvider.enable &&
			!pcfg.alignmentProvider.parametersFilepath.empty()) {
		AlignmentParametersProvider::Config provCfg;
		provCfg.filePath = pcfg.alignmentProvider.parametersFilepath;
		provCfg.treeName = pcfg.alignmentProvider.treeName;

		AlignmentParametersProvider provider(provCfg);
		auto providerStore = provider.getAlignmentStore();

		// merge provider store into aStore (overwrite existing entries)
		for (const auto& entry : *providerStore) {
			(*aStore)[entry.first] = entry.second;
		}
	}

	// Make the store visible to AlignmentAlgorithm’s updater
	g_alignmentStore = aStore;
	// Alignment context
  AlignmentContext alignCtx(aStore);

  // Magnetic Field
  auto field = ApollonGeometry::buildMagField(gctx);

	//track fitting services setup
  using ActionList = KFTrackFittingAlgorithm::ActionList;
	using AbortList  = KFTrackFittingAlgorithm::AbortList;
	using Propagator = FittingServices::Propagator;
	using Trajectory = FittingServices::Trajectory;
	using Stepper    = FittingServices::Stepper;
  using Navigator  = FittingServices::Navigator;

	// Surface accessor for SimpleSourceLink
	SimpleSourceLink::SurfaceAccessor surfaceAccessor(detector.get());

	// Common updater / smoother
	Acts::GainMatrixUpdater  commonUpdater;
	Acts::GainMatrixSmoother commonSmoother;

	// Base extensions
	Acts::KalmanFitterExtensions<Trajectory> baseExtensions;
	baseExtensions.calibrator
			.connect<&simpleSourceLinkCalibrator<Trajectory>>();
	baseExtensions.updater
			.connect<&Acts::GainMatrixUpdater::operator()<Trajectory>>(
					&commonUpdater);
	baseExtensions.smoother
			.connect<&Acts::GainMatrixSmoother::operator()<Trajectory>>(
					&commonSmoother);
	baseExtensions.surfaceAccessor
			.connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
					&surfaceAccessor);

	// Reference surface for sampling the track
	double halfX = std::numeric_limits<double>::max();
	double halfY = std::numeric_limits<double>::max();

	Acts::RotationMatrix3 rotX =
			Acts::AngleAxis3(goInst.toWorldAngleX, Acts::Vector3::UnitX())
					.toRotationMatrix();
	Acts::RotationMatrix3 rotY =
			Acts::AngleAxis3(goInst.toWorldAngleY, Acts::Vector3::UnitY())
					.toRotationMatrix();
	Acts::RotationMatrix3 rotZ =
			Acts::AngleAxis3(goInst.toWorldAngleZ, Acts::Vector3::UnitZ())
					.toRotationMatrix();

	Acts::Transform3 refTrf = Acts::Transform3::Identity();
	refTrf.rotate(rotX);
	refTrf.rotate(rotY);
	refTrf.rotate(rotZ);

	auto refSurface =
			Acts::Surface::makeShared<Acts::PlaneSurface>(
					refTrf,
					std::make_shared<Acts::RectangleBounds>(halfX, halfY));
	Acts::GeometryIdentifier refGeoId;
	refGeoId.setExtra(1);
	refSurface->assignGeometryId(std::move(refGeoId));
  
	// Navigator, Stepper, Propagator, KalmanFitter
	Navigator::Config navCfg;
	navCfg.detector = detector.get();
	navCfg.resolvePassive = false;
	navCfg.resolveMaterial = true;
	navCfg.resolveSensitive = true;

	auto navLogger = Acts::getDefaultLogger(
			"DetectorNavigator", globalLogLevel);
	Navigator kfNavigator(navCfg, std::move(navLogger));

	Stepper kfStepper(std::move(field));

	auto propLogger = Acts::getDefaultLogger(
			"Propagator", globalLogLevel);
	Propagator kfPropagator(
			std::move(kfStepper),
			std::move(kfNavigator),
			std::move(propLogger));
	
	auto kfLogger = Acts::getDefaultLogger(
			"DetectorKalmanFilter", globalLogLevel);
	auto kalmanFitter =
			std::make_shared<
					Acts::KalmanFitter<Propagator, Trajectory>>(
							kfPropagator,
							std::move(kfLogger));

	Acts::PropagatorPlainOptions pOptions(gctx, mctx);

	// Base KF options carry gctx/mctx/cctx and default flags
	Acts::KalmanFitterOptions<Trajectory> baseKfOptions(
			gctx,
			mctx,
			std::cref(cctx),
			baseExtensions,
			pOptions,
			refSurface.get(), /* reference surface */
			true,    /* mScattering */
			true,    /* eLoss */
			false,   /* reverseFiltering */
			1.0      /* rScaling */
	);

	// Set fitting services singleton (to be used by fitting algorithms)
	auto& svc = FittingServices::instance();
	svc.detector        = detector;
	svc.propagator    = std::make_shared<Propagator>(std::move(kfPropagator));
	svc.kalmanFitter  = kalmanFitter;
	svc.baseKfOptions = baseKfOptions;
	svc.baseExtensions = baseExtensions;
	svc.surfaceAccessor.emplace(surfaceAccessor);
	svc.referenceSurface = refSurface;

  // --- FastSim helpers (digitizer + generators) ---
  SimpleDigitizer::Config digitizerCfg;
  digitizerCfg.resolution = {5_um, 5_um};
  svc.simDigitizer = std::make_shared<SimpleDigitizer>(digitizerCfg);

  Acts::Vector3 vertexMean = Acts::Vector3::Zero();
  Acts::SquareMatrix3 vertexCov =
      Acts::SquareMatrix3::Identity() * 0_um;
  svc.simVertexGenerator =
      std::make_shared<GaussianVertexGenerator>(vertexMean, vertexCov);

  svc.simMomentumGenerator = std::make_shared<SphericalMomentumGenerator>();
  svc.simMomentumGenerator->pRange     = std::make_pair(0.1_GeV, 0.7_GeV);
  svc.simMomentumGenerator->thetaRange =
      std::make_pair(M_PI_2 - 5e-3, M_PI_2 + 5e-3);
  svc.simMomentumGenerator->phiRange   =
      std::make_pair(-5e-3, 5e-3);

  svc.simDetSurfaces.clear();
  for (const auto* vol : detector->volumes()) {
    for (const auto* surf : vol->surfaces()) {
      if (surf->geometryId().sensitive()) {
        svc.simDetSurfaces.push_back(surf);
      }
    }
  }
  

  // Sequencer
  Sequencer::Config seqCfg;
  seqCfg.logLevel = globalLogLevel;
  // seqCfg.events   = -1;
  seqCfg.numThreads = 1;
  Sequencer sequencer(seqCfg);

	//geometry context decorator
  // Only decorate geometry if we actually have a misalignment pattern
  if (pcfg.detector.storeMode != "none") {
    sequencer.addContextDecorator(
        std::make_shared<GeometryContextDecorator>(aStore));
  }

	//--------------- build pipeline components ----------------
  // PREPROCESSING
  if (pcfg.preprocessing.enable) {
    pcfg.preprocessing.outputFile = runRoot + "/preprocessed.root";
    // Build PreprocessingConfig
    TrackingPipeline::Preprocessing::PreprocessingConfig preCfg;
    preCfg.enable = true;  // not used inside runPreprocessing, but harmless
    preCfg.inputDirs = pcfg.preprocessing.inputDirs;
    preCfg.inputTreeName = pcfg.preprocessing.inputTreeName;
    preCfg.inputBranchName = pcfg.preprocessing.inputBranchName;
    preCfg.outputFile = pcfg.preprocessing.outputFile;
    preCfg.outputTreeName = pcfg.preprocessing.outputTreeName;
    preCfg.skipEntries = pcfg.preprocessing.skipEntries;
    std::cout << "Running preprocessing on " << preCfg.inputDirs.size()
              << " input directories, writing to '" << preCfg.outputFile << "'\n";
    TrackingPipeline::Preprocessing::runPreprocessing(preCfg);

    try {
      auto& readersection = toml::find(pcfg.root, pcfg.dataReader.type);
      // overwrite / create filepaths array
      readersection.as_table()["filePaths"] =
          toml::array{ pcfg.preprocessing.outputFile };
    } catch (const std::exception& e) {
      std::cerr << "Failed to override reader's filepaths for preprocessing: "
                << e.what() << "\n";
      return 1;
    }
    //debug:
    auto& readerSection = toml::find(pcfg.root, pcfg.dataReader.type);
    auto filePaths = toml::find<std::vector<std::string>>(readerSection, "filePaths");
    std::cout << "DEBUG: reader filePaths after preprocessing:\n";
    for (const auto& p : filePaths) {
      std::cout << "  " << p << "\n";
    }
  }

  // Reader
  try {
    std::cout << "Building reader of type: " << pcfg.dataReader.type << "\n";
    auto reader = TrackingPipeline::buildReader(
        pcfg.dataReader.type,
        pcfg.root,
        surfaceMap,
        globalLogLevel);

    sequencer.addReader(reader);
  } catch (const std::exception& e) {
    std::cerr << "Failed to build reader: " << e.what() << "\n";
    return 1;
  }

  // Seeding algorithms
  try {
    if (pcfg.seeding.enable && !pcfg.seeding.types.empty()) {
      std::cout << "Building " << pcfg.seeding.types.size() << "\n";

      auto seedingAlgos = TrackingPipeline::buildAlgorithms(
          pcfg.seeding.types,
          pcfg.root,
          globalLogLevel);

      for (auto& algo : seedingAlgos) {
        sequencer.addAlgorithm(algo);
      }
    } else {
      std::cout << "Seeding disabled or no types configured\n";
    }
  } catch (const std::exception& e) {
    std::cerr << "Failed to build seeding algorithms: " << e.what() << "\n";
    return 1;
  }

  // Track fitting algorithms
  try {
    if (pcfg.trackFitting.enable && !pcfg.trackFitting.types.empty()) {
      std::cout << "Building " << pcfg.trackFitting.types.size()
                << " track fitting algorithms from " << "\n";

      auto fittingAlgos = TrackingPipeline::buildAlgorithms(
          pcfg.trackFitting.types,
          pcfg.root,
          globalLogLevel);

      for (auto& algo : fittingAlgos) {
        sequencer.addAlgorithm(algo);
      }
    } else {
      std::cout << "Track fitting disabled or no types configured\n";
    }
  } catch (const std::exception& e) {
    std::cerr << "Failed to build track fitting algorithms: " << e.what() << "\n";
    return 1;
  }

  // Track cleaning algorithms
  try {
    if (pcfg.trackCleaning.enable && !pcfg.trackCleaning.types.empty()) {
      std::cout << "Building " << pcfg.trackCleaning.types.size()
                << " track cleaning algorithms\n";

      auto cleaningAlgos = TrackingPipeline::buildAlgorithms(
          pcfg.trackCleaning.types,
          pcfg.root,
          globalLogLevel);

      for (auto& algo : cleaningAlgos) {
        sequencer.addAlgorithm(algo);
      }
    } else {
      std::cout << "Track cleaning disabled or no types configured\n";
    }
  } catch (const std::exception& e) {
    std::cerr << "Failed to build track cleaning algorithms: " << e.what() << "\n";
    return 1;
  }

  // Measurement embedding (FastSim) algorithms
  try {
    if (pcfg.measEmbedding.enable && !pcfg.measEmbedding.types.empty()) {
      std::cout << "Building " << pcfg.measEmbedding.types.size()
                << " measurement embedding algorithms\n";

      auto measAlgos = TrackingPipeline::buildAlgorithms(
          pcfg.measEmbedding.types,
          pcfg.root,
          globalLogLevel);

      for (auto& algo : measAlgos) {
        sequencer.addAlgorithm(algo);
      }
    } else {
      std::cout << "Measurement embedding disabled or no types configured\n";
    }
  } catch (const std::exception& e) {
    std::cerr << "Failed to build measurement embedding algorithms: "
              << e.what() << "\n";
    return 1;
  }

  // Writers
  try {
    if (!pcfg.dataWriter.types.empty()) {
      std::cout << "Building " << pcfg.dataWriter.types.size() << "\n";

      auto writers = TrackingPipeline::buildWriters(
          pcfg.dataWriter.types,
          pcfg.root,
          globalLogLevel,
          runRoot);

      for (auto& w : writers) {
        sequencer.addWriter(w);
      }
    } else {
      std::cout << "No writers configured (DATA_WRITER.types is empty)\n";
    }
  } catch (const std::exception& e) {
    std::cerr << "Failed to build writers: " << e.what() << "\n";
    return 1;
  }

  // Run
  std::cout << "Starting sequencer run...\n";
  try {
    sequencer.run();
  } catch (const std::exception& e) {
    std::cerr << "Sequencer runtime error: " << e.what() << "\n";
    return 1;
  }

  std::cout << "Run finished successfully.\n";
  return 0;
}

} // namespace TrackingPipeline

