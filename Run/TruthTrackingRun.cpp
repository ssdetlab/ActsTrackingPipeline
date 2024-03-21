#include "ActsLUXEPipeline/LUXEROOTDataReader.hpp"
#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXEIdealSeeder.hpp"
#include "ActsLUXEPipeline/LUXETrackFitter.hpp"
#include "ActsLUXEPipeline/Sequencer.hpp"
#include "ActsLUXEPipeline/EventVisualizer.hpp" 
#include "ActsLUXEPipeline/LUXEMagneticField.hpp"

#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"

#include <filesystem>

#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

int main() {
    // Set the log level
    Acts::Logging::Level logLevel = Acts::Logging::VERBOSE;

    // Dummy context and options
    Acts::GeometryContext gctx;
    Acts::MagneticFieldContext mctx;
    Acts::CalibrationContext cctx;
    LUXEGeometry::GeometryOptions gOpt;

    // Set the path to the gdml file
    // and the names of the volumes to be converted
    std::string gdmlPath = 
        "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_gdmls/lxgeomdump_stave_positron.gdml";
    std::vector<std::string> names{"OPPPSensitive"};

    // Build the LUXE detector
    auto trackerBP = 
        LUXEGeometry::makeBlueprintLUXE(gdmlPath, names, gOpt);
    auto detector =
        LUXEGeometry::buildLUXEDetector(std::move(trackerBP), gctx, gOpt);

    for (auto& vol : detector->rootVolumes()) {
        for (auto surf : vol->surfaces()) {
            std::cout << "Surface " << surf->geometryId() << " " << surf->center(gctx).transpose() << std::endl;
        }
    }

    // Setup the sequencer
    Sequencer::Config seqCfg;
    seqCfg.events = 10;
    seqCfg.numThreads = -1;
    Sequencer sequencer(seqCfg);

    // Add the sim data reader
    LUXEROOTReader::LUXEROOTSimDataReader::Config readerCfg 
        = LUXEROOTReader::defaultSimConfig();
    readerCfg.dataCollection = "SimMeasurements";
    std::string pathToDir = 
        "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_dataInRootFormat/SignalNextTrial_e1gpc_10.0_1";

    // Get the paths to the files in the directory
    for (const auto & entry : std::filesystem::directory_iterator(pathToDir)) {
        std::string pathToFile = entry.path();
        readerCfg.filePaths.push_back(pathToFile);
    }

    // The events are not sorted in the directory
    // but we need to process them in order
    std::sort(readerCfg.filePaths.begin(), readerCfg.filePaths.end(),
        [] (const std::string& a, const std::string& b) {
            std::size_t idxRootA = a.find_last_of('.');
            std::size_t idxEventA = a.find_last_of('t', idxRootA);
            std::string eventSubstrA = a.substr(idxEventA + 1, idxRootA - idxEventA);
            
            std::size_t idxRootB = b.find_last_of('.');
            std::size_t idxEventB = b.find_last_of('t', idxRootB);
            std::string eventSubstrB = b.substr(idxEventB + 1, idxRootB - idxEventB);

            return std::stoul(eventSubstrA) < std::stoul(eventSubstrB);
        }
    );

    // Process only the first event
    readerCfg.filePaths = std::vector<std::string>(
        readerCfg.filePaths.begin(), readerCfg.filePaths.begin() + 72); 

    // Add the reader to the sequencer
    sequencer.addReader(
        std::make_shared<LUXEROOTReader::LUXEROOTSimDataReader>(readerCfg, logLevel));

    // Add the ideal seeder to the sequencer
    IdealSeeder::Config seederCfg{
        .inputCollection = "SimMeasurements",
        .outputCollection = "SimSeeds",
        .minHits = 4,
        .maxHits = 4,
        .gOpt = gOpt
    };
    sequencer.addAlgorithm(
        std::make_shared<IdealSeeder>(seederCfg, logLevel));

    // // Add the event visualizer to the sequencer
    // EventVisualizer::Config visCfg{
        // .inputCollection = "SimSeeds",
        // .outputPath = "hits.obj",
        // .nTracks = 1,
        // .detector = detector.get(),
        // .visualizeVolumes = false,
        // .visualizeHits = true,
        // .visualizeTracks = false
    // };

    // sequencer.addAlgorithm(
        // std::make_shared<EventVisualizer>(visCfg, logLevel));

    Acts::GainMatrixUpdater kfUpdater;
    Acts::GainMatrixSmoother kfSmoother;

    SimpleSourceLink::SurfaceAccessor surfaceAccessor{*detector};

    // Initialize track fitter options
    Acts::KalmanFitterExtensions<
        Acts::VectorMultiTrajectory> extensions;
    // Add calibrator
    extensions.calibrator.connect<&simpleSourceLinkCalibrator<
        Acts::VectorMultiTrajectory>>();
    // Add the updater
    extensions.updater.connect<&Acts::GainMatrixUpdater::operator()<
            Acts::VectorMultiTrajectory>>(&kfUpdater);
    // Add the smoother
    extensions.smoother.connect<&Acts::GainMatrixSmoother::operator()<
        Acts::VectorMultiTrajectory>>(&kfSmoother);
    // Add the surface accessor
    extensions.surfaceAccessor.connect<
        &SimpleSourceLink::SurfaceAccessor::operator()>(
        &surfaceAccessor);

    using ActionList = Acts::ActionList<>;
    using AbortList = Acts::AbortList<Acts::EndOfWorldReached>;

    using Propagator = Acts::Propagator<Acts::EigenStepper<>, 
        Acts::Experimental::DetectorNavigator>;

    using KF = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;

    auto propOptions = 
        Acts::PropagatorOptions<ActionList,AbortList>(gctx, mctx);

    propOptions.maxSteps = 20;

    auto options = Acts::KalmanFitterOptions(gctx, mctx, cctx, extensions,
        propOptions);

    using namespace Acts::UnitLiterals;
    Acts::Experimental::DetectorNavigator::Config cfg;
    cfg.detector = detector.get();
    cfg.resolvePassive = false;
    cfg.resolveMaterial = true;
    cfg.resolveSensitive = true;
    Acts::Experimental::DetectorNavigator navigator(
        cfg, Acts::getDefaultLogger("DetectorNavigator", Acts::Logging::VERBOSE));

    // const std::vector<std::pair<
        // Acts::ActsScalar,Acts::ActsScalar>> MagneticFieldBounds =
            // {std::make_pair(-255_mm,255_mm),
            // std::make_pair(-255_mm,255_mm),
            // std::make_pair(-255_mm,255_mm)};

    // auto transformPos = [&](const Acts::Vector3& pos) {
        // for (int i=0;i<3;i++) {
            // if (pos[i] < MagneticFieldBounds[i].first ||
                // pos[i] > MagneticFieldBounds[i].second) {
                // return Acts::Vector3{0,0,0};
            // }
        // }
        // return pos;
    // };

    // auto transformBField = [](const Acts::Vector3& field, const Acts::Vector3&) {
        // return field;
    // };

    // const std::vector<unsigned int> bins{5u, 5u, 5u};

    // auto bFieldVal = LUXEMagneticField::SimpleDipole3(); 
    // bFieldVal.m_magnitude = 1.2_T;
    // auto BField = 
        // LUXEMagneticField::buildBFieldMapEq3(
            // bFieldVal, 
            // transformPos, transformBField, bins);

    // auto BFieldPtr = std::make_shared<decltype(BField)>(BField);

    auto field =
        std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, 0.0));

    Acts::EigenStepper<> stepper(std::move(field));
    auto propagator = Propagator(
        std::move(stepper), std::move(navigator),
        Acts::getDefaultLogger("Propagator", Acts::Logging::VERBOSE));

    auto dkfLogger = Acts::getDefaultLogger("DetectorKalmanFilter", Acts::Logging::VERBOSE);
    const auto fitter = KF(propagator, std::move(dkfLogger));

    // Add the track fitting algorithm to the sequencer
    TrackFitter<Propagator>::Config fitterCfg{
        .inputCollection = "SimSeeds",
        .outputCollection = "SimTracks",
        .fitter = fitter,
        .kfOptions = options};

    sequencer.addAlgorithm(
        std::make_shared<TrackFitter<Propagator>>(fitterCfg, logLevel));

    // Run all configured algorithms and return the appropriate status.
    return sequencer.run();
}