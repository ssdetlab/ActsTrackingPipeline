#include "ActsLUXEPipeline/E320ROOTDataReader.hpp"
#include "ActsLUXEPipeline/E320Geometry.hpp"
#include "ActsLUXEPipeline/TrackFittingAlgorithm.hpp"
#include "ActsLUXEPipeline/QuadrupoleMagField.hpp"
#include "ActsLUXEPipeline/DipoleMagField.hpp"
#include "ActsLUXEPipeline/CompositeMagField.hpp"
#include "ActsLUXEPipeline/Sequencer.hpp"
#include "ActsLUXEPipeline/ROOTFittedTrackWriter.hpp"
#include "ActsLUXEPipeline/E320SourceLinkGrid.hpp"
#include "ActsLUXEPipeline/CsvLookupTableProvider.hpp"
#include "ActsLUXEPipeline/ForwardOrderedIntersectionFinder.hpp"
#include "ActsLUXEPipeline/E320PathWidthProvider.hpp"
#include "ActsLUXEPipeline/PathSeedingAlgorithm.hpp"
#include "ActsLUXEPipeline/TryAllTrackFindingAlgorithm.hpp"
#include "ActsLUXEPipeline/CKFTrackFindingAlgorithm.hpp"
#include "ActsLUXEPipeline/Generators.hpp"
#include "ActsLUXEPipeline/NoiseEmbeddingAlgorithm.hpp"
#include "ActsLUXEPipeline/PhoenixTrackWriter.hpp"

#include "Acts/Seeding/PathSeeder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"

#include <filesystem>

using Grid = E320TrackFinding::E320SourceLinkGridConstructor::GridType;

using ActionList = Acts::ActionList<>;
using AbortList = Acts::AbortList<Acts::EndOfWorldReached>;

using Propagator = Acts::Propagator<
    Acts::EigenStepper<>, 
    Acts::Experimental::DetectorNavigator>;
using PropagatorOptions =
    typename Propagator::template Options<ActionList, AbortList>;

using Trajectory = Acts::VectorMultiTrajectory;
using KFTrackContainer = Acts::VectorTrackContainer;
using KF = Acts::KalmanFitter<Propagator, Trajectory>;

using CKFTrackContainer = Acts::TrackContainer<
    Acts::VectorTrackContainer,
    Acts::VectorMultiTrajectory,
    Acts::detail::ValueHolder>;

using TrackStateContainerBackend =
    typename CKFTrackContainer::TrackStateContainerBackend;

using namespace Acts::UnitLiterals;

int main() {
    // Set the log level
    Acts::Logging::Level logLevel = Acts::Logging::FATAL;

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
        "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/E320Pipeline_gdmls/ettgeom_magnet_pdc_tracker.gdml";
    std::vector<std::string> names{"OPPPSensitive", "DetChamberWindow"};

    std::string materialPath = "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/E320Pipeline_material/uniform/material.json";

    // Build the detector
    auto trackerBP = 
        E320Geometry::makeBlueprintE320(gdmlPath, names, gOpt);
    auto detector =
        E320Geometry::buildE320Detector(std::move(trackerBP), gctx, gOpt, materialPath, {});

    for (auto& v : detector->volumes()) {
        std::cout << v->name() << std::endl;
        for (auto& s : v->surfaces()) {
            std::cout << s->center(gctx).transpose() << "   " << s->normal(gctx, s->center(gctx), Acts::Vector3(1,0,0)).transpose() << std::endl;
        }
    }

    // --------------------------------------------------------------
    // The magnetic field setup

    // Extent in already rotated frame
    Acts::Extent quad1Extent;
    quad1Extent.set(
        Acts::BinningValue::binX, 
        gOpt.quad1Translation[0] - gOpt.quad1Bounds[0],
        gOpt.quad1Translation[0] + gOpt.quad1Bounds[0]);
    quad1Extent.set(
        Acts::BinningValue::binZ,
        gOpt.quad1Translation[1] - gOpt.quad1Bounds[1],
        gOpt.quad1Translation[1] + gOpt.quad1Bounds[1]);
    quad1Extent.set(
        Acts::BinningValue::binY,
        gOpt.quad1Translation[2] - gOpt.quad1Bounds[2],
        gOpt.quad1Translation[2] + gOpt.quad1Bounds[2]);

    Acts::Extent quad2Extent;
    quad2Extent.set(
        Acts::BinningValue::binX, 
        gOpt.quad2Translation[0] - gOpt.quad2Bounds[0],
        gOpt.quad2Translation[0] + gOpt.quad2Bounds[0]);
    quad2Extent.set(
        Acts::BinningValue::binZ,
        gOpt.quad2Translation[1] - gOpt.quad2Bounds[1],
        gOpt.quad2Translation[1] + gOpt.quad2Bounds[1]);
    quad2Extent.set(
        Acts::BinningValue::binY,
        gOpt.quad2Translation[2] - gOpt.quad2Bounds[2],
        gOpt.quad2Translation[2] + gOpt.quad2Bounds[2]);

    Acts::Extent quad3Extent;
    quad3Extent.set(
        Acts::BinningValue::binX, 
        gOpt.quad3Translation[0] - gOpt.quad3Bounds[0],
        gOpt.quad3Translation[0] + gOpt.quad3Bounds[0]);
    quad3Extent.set(
        Acts::BinningValue::binZ,
        gOpt.quad3Translation[1] - gOpt.quad3Bounds[1],
        gOpt.quad3Translation[1] + gOpt.quad3Bounds[1]);
    quad3Extent.set(
        Acts::BinningValue::binY,
        gOpt.quad3Translation[2] - gOpt.quad3Bounds[2],
        gOpt.quad3Translation[2] + gOpt.quad3Bounds[2]);

    Acts::Extent dipoleExtent;
    dipoleExtent.set(
        Acts::BinningValue::binX, 
        gOpt.dipoleTranslation.x() - gOpt.dipoleBounds[0],
        gOpt.dipoleTranslation.x() + gOpt.dipoleBounds[0]);
    dipoleExtent.set(
        Acts::BinningValue::binZ,
        gOpt.dipoleTranslation.y() - gOpt.dipoleBounds[1],
        gOpt.dipoleTranslation.y() + gOpt.dipoleBounds[1]);
    dipoleExtent.set(
        Acts::BinningValue::binY,
        gOpt.dipoleTranslation.z() - gOpt.dipoleBounds[2],
        gOpt.dipoleTranslation.z() + gOpt.dipoleBounds[2]);

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
    
    Acts::ActsScalar dipoleB = 0.31_T;
    DipoleMagField dipoleField(
        gOpt.dipoleParams, 
        dipoleB,
        gOpt.actsToWorldRotation, 
        gOpt.actsToWorldRotation.inverse() * gOpt.dipoleTranslation);
    
    CompositeMagField::FieldComponents fieldComponents = {
        {quad1Extent, &quad1Field},
        {quad2Extent, &quad2Field},
        {quad3Extent, &quad3Field},
        {dipoleExtent, &dipoleField}
    };

    auto field = std::make_shared<CompositeMagField>(fieldComponents);

    // --------------------------------------------------------------
    // Event reading 

    // Setup the sequencer
    Sequencer::Config seqCfg;
    // seqCfg.events = 2;
    seqCfg.numThreads = -1;
    seqCfg.trackFpes = false;
    Sequencer sequencer(seqCfg);

    // Add the sim data reader
    E320ROOTReader::E320ROOTSimDataReader::Config readerCfg = 
        E320ROOTReader::defaultSimConfig();
    readerCfg.outputSourceLinks = "Measurements";
    readerCfg.outputSimClusters = "SimClusters";
    std::string pathToDir = 
        "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/E320Pipeline_dataInRootFormat/Signal_E320lp_10.0_12BX_All/Signal_E320lp_10.0_12BX_10us_integration_time";

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

    // Add the reader to the sequencer
    sequencer.addReader(
        std::make_shared<E320ROOTReader::E320ROOTSimDataReader>(readerCfg, logLevel));

    // --------------------------------------------------------------
    // The path seeding setup
    SimpleSourceLink::SurfaceAccessor surfaceAccessor{detector.get()};

    auto pathSeederCfg = Acts::Experimental::PathSeeder::Config();

    // Estimator of the IP and first hit
    // parameters of the track
    CsvLookupTableProvider::Config trackLookupCfg;

    trackLookupCfg.filePath = 
        "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/E320Pipeline_lookups/RangedUniform_05_45_Stationary_000_200k_1000x100_MaterialOn_lookup.csv";

    CsvLookupTableProvider trackLookup(trackLookupCfg);
    pathSeederCfg.trackEstimator.connect<
        &CsvLookupTableProvider::operator()>(
        &trackLookup);

    // Transforms the source links to global coordinates
    SimpleSourceLinkCoordinateCalibrator sourceLinkCalibrator;
    sourceLinkCalibrator.m_surfaceAccessor.connect<
        &SimpleSourceLink::SurfaceAccessor::operator()>(
        &surfaceAccessor);
    pathSeederCfg.sourceLinkCalibrator.connect<
        &SimpleSourceLinkCoordinateCalibrator::operator()>(
        &sourceLinkCalibrator);

    // Intersection finder
    ForwardOrderedIntersectionFinder intersectionFinder;

    // Combine layers into surfaces for the intersection finder
    std::vector<std::shared_ptr<Acts::Surface>> surfacePtrs;
    for (int i = 0; i < gOpt.staveZ.size(); i++) {
        double halfX = gOpt.chipSizeX/2;
        double halfY = ((gOpt.chipY.at(8) + gOpt.chipSizeX/2) - 
            (gOpt.chipY.at(0) - gOpt.chipSizeX/2))/2;

        double centerY = ((-gOpt.chipY.at(8) - gOpt.chipSizeY/2) + 
            (-gOpt.chipY.at(0) + gOpt.chipSizeY/2))/2;

        Acts::Transform3 transform(
            Acts::Translation3(Acts::Vector3(gOpt.chipX, gOpt.staveZ.at(i), centerY)) *
            gOpt.actsToWorldRotation.inverse());

        auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
            transform,
            std::make_shared<Acts::RectangleBounds>(
                halfX, halfY));

        Acts::GeometryIdentifier geoId;
        geoId.setSensitive(i + 1);
        surface->assignGeometryId(std::move(geoId));
        surfacePtrs.push_back(surface);

        if (i == 0) {
            pathSeederCfg.firstLayerIds.push_back(geoId);
        }
    }
    std::vector<const Acts::Surface*> surfaces;
    for (auto surface : surfacePtrs) {
        surfaces.push_back(surface.get());
    }

    intersectionFinder.m_surfaces = std::move(surfaces);
    intersectionFinder.m_tol = 
        (gOpt.chipY.at(1) - gOpt.chipSizeY/2) - (gOpt.chipY.at(0) + gOpt.chipSizeY/2) + 1_mm;

    pathSeederCfg.intersectionFinder.connect<
        &ForwardOrderedIntersectionFinder::operator()>(&intersectionFinder);

    // Path width provider
    std::map<std::int32_t, std::pair<
        Acts::ActsScalar,Acts::ActsScalar>> 
            pathWidths = {
                {0, {100_um, 100_um}},
                {1, {200_um, 200_um}},
                {2, {250_um, 250_um}},
                {3, {300_um, 300_um}},
    };

    E320TrackFinding::E320PathWidthProvider pathWidthProvider(
        gOpt,
        pathWidths);

    pathSeederCfg.pathWidthProvider.connect<
        &E320TrackFinding::E320PathWidthProvider::operator()>(
        &pathWidthProvider);

    pathSeederCfg.orientation = Acts::BinningValue::binY;

    // Grid to bin the source links
    E320TrackFinding::E320SourceLinkGridConstructor::Config gridConstructorCfg{
        .gOpt = gOpt,
        .bins = std::make_pair(20, 1000),
    };
    gridConstructorCfg.surfaceAccessor.connect<
        &SimpleSourceLink::SurfaceAccessor::operator()>(
        &surfaceAccessor);

    auto gridConstructor = std::make_shared<E320TrackFinding::E320SourceLinkGridConstructor>(gridConstructorCfg); 

    // Create the path seeder algorithm
    auto seedingAlgoCfg = PathSeedingAlgorithm::Config();
    seedingAlgoCfg.seeder = std::make_shared<Acts::Experimental::PathSeeder>(pathSeederCfg);
    seedingAlgoCfg.sourceLinkGridConstructor = gridConstructor;
    seedingAlgoCfg.inputSourceLinks = "Measurements";
    seedingAlgoCfg.outputSeeds = "PathSeeds";

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
    auto ckfPropagator = Propagator(
        std::move(ckfStepper), std::move(ckfNavigator),
        Acts::getDefaultLogger("Propagator", logLevel));

    Acts::CombinatorialKalmanFilter<Propagator,CKFTrackContainer> ckf(
        ckfPropagator, Acts::getDefaultLogger("CombinatorialKalmanFilter", logLevel));

    // Configuration for the measurement selector
    std::vector<std::pair<
        Acts::GeometryIdentifier, Acts::MeasurementSelectorCuts>> cuts;
    for (auto& vol : detector->volumes()) {
        for (auto& surf : vol->surfaces()) {
            if (vol->name() == "layer0") {
                cuts.push_back(
                    {
                        surf->geometryId(), 
                        {
                            {}, 
                            {std::numeric_limits<Acts::ActsScalar>::max()}, 
                            {1000u}
                        }
                    });
            }
            else {
                cuts.push_back(
                    {
                        surf->geometryId(), 
                        {
                            {}, 
                            {1}, 
                            {1u}
                        }
                    });
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
        &Acts::GainMatrixUpdater::operator()<TrackStateContainerBackend>>(&ckfUpdater);
    ckfExtensions.measurementSelector.template connect<
        &Acts::MeasurementSelector::select<TrackStateContainerBackend>>(
        &measSel);

    CKFTrackFindingAlgorithm<Propagator, CKFTrackContainer>::Config trackFindingCfg{
        .ckf = ckf,
    };
    trackFindingCfg.extensions = ckfExtensions;
    trackFindingCfg.inputSeeds = "PathSeeds";
    trackFindingCfg.outputTrackCandidates = "TrackCandidates";
    trackFindingCfg.minSourceLinks = 4;
    trackFindingCfg.maxSourceLinks = 4;

    auto trackFindingAlgorithm = 
        std::make_shared<CKFTrackFindingAlgorithm<Propagator, CKFTrackContainer>>(trackFindingCfg, logLevel);
    sequencer.addAlgorithm(trackFindingAlgorithm);

    // --------------------------------------------------------------
    // Track fitting

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
        PropagatorOptions(gctx, mctx);

    propOptions.maxSteps = 300;

    auto options = Acts::KalmanFitterOptions(gctx, mctx, cctx, extensions,
        propOptions);

    // Reference surface for sampling the track at the IP
    double halfX = std::numeric_limits<double>::max();
    double halfY = std::numeric_limits<double>::max();

    Acts::Transform3 transform(
        Acts::Translation3(Acts::Vector3(0, 0, 0)) *
        gOpt.actsToWorldRotation.inverse());

    auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
        transform,
        std::make_shared<Acts::RectangleBounds>(
            halfX, halfY));

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
    auto kfPropagator = Propagator(
        std::move(kfStepper), std::move(kfNavigator),
        Acts::getDefaultLogger("Propagator", logLevel));

    const auto fitter = 
        KF(kfPropagator, 
            Acts::getDefaultLogger("DetectorKalmanFilter", logLevel));

    // Add the track fitting algorithm to the sequencer
    TrackFittingAlgorithm<
        Propagator, 
        Trajectory, 
        KFTrackContainer>::Config fitterCfg{
            .inputCollection = "TrackCandidates",
            .outputCollection = "Tracks",
            .fitter = fitter,
            .kfOptions = options};

    sequencer.addAlgorithm(
        std::make_shared<
            TrackFittingAlgorithm<
            Propagator, 
            Trajectory, 
            KFTrackContainer>>(fitterCfg, logLevel));

    // --------------------------------------------------------------
    // Event write out

    auto trackWriterCfg = ROOTFittedTrackWriter::Config();
    trackWriterCfg.surfaceAccessor.connect<
        &SimpleSourceLink::SurfaceAccessor::operator()>(
            &surfaceAccessor);

    trackWriterCfg.inputKFTracks = "Tracks";
    trackWriterCfg.inputTruthClusters = "SimClusters";
    trackWriterCfg.treeName = "fitted-tracks";
    trackWriterCfg.filePath = "fitted-tracks-sig-2k.root";

    sequencer.addWriter(
        std::make_shared<ROOTFittedTrackWriter>(trackWriterCfg, logLevel));

    // // --------------------------------------------------------------
    // // Phoenix write out
    // PhoenixTrackWriter::Config phoenixWriterCfg;
    // phoenixWriterCfg.inputTrackCollection = "Tracks";
    // phoenixWriterCfg.fileName = "test-tracks";

    // sequencer.addWriter(
        // std::make_shared<PhoenixTrackWriter>(phoenixWriterCfg, logLevel));

    return sequencer.run();
}
