#include "ActsLUXEPipeline/LUXEROOTDataReader.hpp"
#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/IdealSeeder.hpp"
#include "ActsLUXEPipeline/TrackFitter.hpp"
#include "ActsLUXEPipeline/ConstantBoundedField.hpp"
#include "ActsLUXEPipeline/Sequencer.hpp"
#include "ActsLUXEPipeline/ROOTFittedTrackWriter.hpp"

#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"

#include <filesystem>

#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

using ActionList = Acts::ActionList<>;
using AbortList = Acts::AbortList<Acts::EndOfWorldReached>;

using Propagator = Acts::Propagator<
    Acts::EigenStepper<>, 
    Acts::Experimental::DetectorNavigator>;

using Trajectory = Acts::VectorMultiTrajectory;
using TrackContainer = Acts::VectorTrackContainer;
using KF = Acts::KalmanFitter<Propagator, Trajectory>;

using namespace Acts::UnitLiterals;

int main() {
    // Set the log level
    Acts::Logging::Level logLevel = Acts::Logging::INFO;

    // Dummy context and options
    Acts::GeometryContext gctx;
    Acts::MagneticFieldContext mctx;
    Acts::CalibrationContext cctx;
    LUXEGeometry::GeometryOptions gOpt;

    // --------------------------------------------------------------
    // LUXE detector setup

    // Set the path to the gdml file
    // and the names of the volumes to be converted
    std::string gdmlPath = 
        "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_gdmls/lxgeomdump_stave_positron.gdml";
    std::vector<std::string> names{"OPPPSensitive"};

    std::string materialPath = "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_build/material.json";

    // Build the LUXE detector
    auto trackerBP = 
        LUXEGeometry::makeBlueprintLUXE(gdmlPath, names, gOpt);
    auto detector =
        LUXEGeometry::buildLUXEDetector(std::move(trackerBP), gctx, gOpt, materialPath, {});

    // --------------------------------------------------------------
    // The magnetic field setup

    // Extent in already rotated frame
    Acts::Extent dipoleExtent;
    dipoleExtent.set(
        Acts::binX, 
        gOpt.dipoleTranslation[0] - gOpt.dipoleBounds[0] + gOpt.constantFieldDelta[0],
        gOpt.dipoleTranslation[0] + gOpt.dipoleBounds[0] - gOpt.constantFieldDelta[0]);
    dipoleExtent.set(
        Acts::binZ,
        gOpt.dipoleTranslation[1] - gOpt.dipoleBounds[1] + gOpt.constantFieldDelta[1],
        gOpt.dipoleTranslation[1] + gOpt.dipoleBounds[1] - gOpt.constantFieldDelta[1]);
    dipoleExtent.set(
        Acts::binY,
        gOpt.dipoleTranslation[2] - gOpt.dipoleBounds[2] + gOpt.constantFieldDelta[2],
        gOpt.dipoleTranslation[2] + gOpt.dipoleBounds[2] - gOpt.constantFieldDelta[2]);

    auto field = std::make_shared<ConstantBoundedField>(
        Acts::Vector3(0., 0., -1.2_T),
        dipoleExtent);

    // --------------------------------------------------------------
    // Event reading 

    // Setup the sequencer
    Sequencer::Config seqCfg;
    seqCfg.events = 10;
    seqCfg.numThreads = 16;
    seqCfg.trackFpes = false;
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
    readerCfg.energyCuts = {1_GeV, 100_GeV};

    // Vertex position extent in the already rotated frame
    Acts::Extent vertexExtent;
    vertexExtent.set(
        Acts::binX, -150_cm, 150_cm);
    vertexExtent.set(
        Acts::binZ, -150_cm, 150_cm);
    vertexExtent.set(
        Acts::binY, -150_cm, 150_cm);

    readerCfg.vertexPosExtent = vertexExtent;

    // Add the reader to the sequencer
    sequencer.addReader(
        std::make_shared<LUXEROOTReader::LUXEROOTSimDataReader>(readerCfg, logLevel));


    // --------------------------------------------------------------
    // Seed finding 

    // Add the ideal seeder to the sequencer
    IdealSeeder::Config seederCfg{
        .inputCollection = "SimMeasurements",
        .outputCollection = "SimSeeds",
        .minHits = 3,
        .maxHits = 10,
    };
    sequencer.addAlgorithm(
        std::make_shared<IdealSeeder>(seederCfg, logLevel));


    // --------------------------------------------------------------
    // Track fitting

    SimpleSourceLink::SurfaceAccessor surfaceAccessor{*detector};

    Acts::GainMatrixUpdater kfUpdater;
    Acts::GainMatrixSmoother kfSmoother;

    // Initialize track fitter options
    Acts::KalmanFitterExtensions<Trajectory> extensions;
    // Add calibrator
    extensions.calibrator.connect<&simpleSourceLinkCalibrator<Trajectory>>();
    // Add the updater
    extensions.updater.connect<
        &Acts::GainMatrixUpdater::operator()<Trajectory>>(&kfUpdater);
    // Add the smoother
    extensions.smoother.connect<
        &Acts::GainMatrixSmoother::operator()<Trajectory>>(&kfSmoother);
    // Add the surface accessor
    extensions.surfaceAccessor.connect<
        &SimpleSourceLink::SurfaceAccessor::operator()>(
        &surfaceAccessor);

    auto propOptions = 
        Acts::PropagatorOptions<ActionList,AbortList>(gctx, mctx);

    propOptions.maxSteps = 200;

    auto options = Acts::KalmanFitterOptions(gctx, mctx, cctx, extensions,
        propOptions);

    Acts::Experimental::DetectorNavigator::Config cfg;
    cfg.detector = detector.get();
    cfg.resolvePassive = false;
    cfg.resolveMaterial = true;
    cfg.resolveSensitive = true;
    Acts::Experimental::DetectorNavigator navigator(
        cfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));

    Acts::EigenStepper<> stepper(std::move(field));
    auto propagator = Propagator(
        std::move(stepper), std::move(navigator),
        Acts::getDefaultLogger("Propagator", logLevel));

    const auto fitter = 
        KF(propagator, 
            Acts::getDefaultLogger("DetectorKalmanFilter", logLevel));

    // Add the track fitting algorithm to the sequencer
    TrackFitter<
        Propagator, 
        Trajectory, 
        TrackContainer>::Config fitterCfg{
            .inputCollection = "SimSeeds",
            .outputCollection = "SimTracks",
            .fitter = fitter,
            .kfOptions = options};

    sequencer.addAlgorithm(
        std::make_shared<
            TrackFitter<
            Propagator, 
            Trajectory, 
            TrackContainer>>(fitterCfg, logLevel));


    // --------------------------------------------------------------
    // Event write out

    auto trackWriterCfg = ROOTFittedTrackWriter::Config{
        "SimTracks",
        "fitted-tracks",
        "fitted-tracks.root",
        3,
        10
    };

    sequencer.addWriter(
        std::make_shared<ROOTFittedTrackWriter>(trackWriterCfg, logLevel));

    // --------------------------------------------------------------
    // Run all configured algorithms and return the appropriate status.

    return sequencer.run();
}
