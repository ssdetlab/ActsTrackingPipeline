#include "ActsLUXEPipeline/E320ROOTDataReader.hpp"
#include "ActsLUXEPipeline/E320Geometry.hpp"
#include "ActsLUXEPipeline/TrackFitter.hpp"
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
#include "ActsLUXEPipeline/IdealSeedingAlgorithm.hpp"

#include "Acts/Seeding/PathSeeder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

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
using TrackContainer = Acts::VectorTrackContainer;
using KF = Acts::KalmanFitter<Propagator, Trajectory>;

using namespace Acts::UnitLiterals;

// --------------------------------------------------------------
// Dummy for comparing the seeders
// Will be removed in the future
class SeedComparer : public IAlgorithm {
    public:
        struct Config {
            std::string inputIdealCollection;
            std::string inputPathCollection;
            Acts::SourceLinkSurfaceAccessor surfaceAccessor;
        };

        SeedComparer(const Config& config, Acts::Logging::Level level) 
        : IAlgorithm("SeedComparer", level), 
        m_cfg(config) {
            m_inputIdealSeeds.initialize(m_cfg.inputIdealCollection);
            m_inputPathSeeds.initialize(m_cfg.inputPathCollection);
        };

        ~SeedComparer() = default;

        ProcessCode execute(const AlgorithmContext& context) const override {
            auto idealSeeds = m_inputIdealSeeds(context);
            auto pathSeeds = m_inputPathSeeds(context);

            if (idealSeeds.empty() || pathSeeds.empty()) {
                return ProcessCode::SUCCESS;
            }

            std::cout << "Ideal seeds size: " << idealSeeds.size() << std::endl;
            std::cout << "Path seeds size: " << pathSeeds.size() << std::endl;

            // Path seeds that contain the pivot
            // from the ideal seeds
            int matchedPivots = 0;
            
            // Path seeds that do not contain
            // the pivot from the ideal seeds
            int lostPivots = 0;

            // Path seeds that arise from 
            // the pivot that is not in the ideal seeds
            int fakePivots = 0;

            // Path seeds that contain 
            // all the measuremets from the ideal seeds
            int matchedSeeds = 0;

            // Path seeds that have one or more
            // measurements lost from the ideal seeds
            int lostSeeds = 0;

            // Normalization
            int idealPivots = idealSeeds.size();

            // Path sees size
            double pathSeedSize = 0;

            E320Geometry::GeometryOptions gOpt;

            for (auto idealSeed : idealSeeds) {
                auto sls = idealSeed.sourceLinks;
                std::sort(sls.begin(), sls.end(),
                    [&](const Acts::SourceLink& a, const Acts::SourceLink& b) {
                        auto sslA = a.get<SimpleSourceLink>();
                        auto sslB = b.get<SimpleSourceLink>();
                        std::size_t idA = static_cast<int>(sslA.geometryId().sensitive()/10 - 1);
                        std::size_t idB = static_cast<int>(sslB.geometryId().sensitive()/10 - 1);
                        return gOpt.staveZ.at(idA) < gOpt.staveZ.at(idB);
                });

                // Select the pivot source link
                auto pivotSl = sls.at(0);

                auto id = static_cast<int>(pivotSl.get<SimpleSourceLink>().geometryId().sensitive()/10 - 1);

                if (id != 0) {
                    idealPivots--;
                    continue;
                }

                // Find the path seed that contains the pivot
                auto pathSeedIt = std::find_if(
                    pathSeeds.begin(), 
                    pathSeeds.end(), 
                    [&](const auto& ps) {
                        // Get the pivot source link
                        // of the ideal seed
                        auto ssl = pivotSl.get<SimpleSourceLink>();
                        
                        // Get the pivot source link
                        // of the path seed
                        Acts::SourceLink sl1 = ps.sourceLinks.at(0);
                        auto ssl1 = sl1.get<SimpleSourceLink>();

                        auto sid = static_cast<int>(ssl1.geometryId().sensitive()/10 - 1);

                        if (sid != 0) {
                            // std::cout << "ALERT " << sid << std::endl;
                        }

                        // Compare the pivot source links
                        return ssl == ssl1;
                });
                // Pivot is not in the ideal seed
                // so the path seed is fake
                if (pathSeedIt == pathSeeds.end()) {
                    // std::cout << "IDEAL PIVOT " << m_cfg.surfaceAccessor(pivotSl)->localToGlobal(context.geoContext, pivotSl.get<SimpleSourceLink>().parameters, Acts::Vector3(0,1,0)).transpose() << std::endl;
                    // std::cout << idealSeed.ipParameters << std::endl;
                    // std::cout << "ENERGY " << idealSeed.ipParameters.absoluteMomentum() << std::endl;
                    lostPivots++;
                    continue;
                }
            }

            // Iterate over the path seeds
            for (auto pathSeed : pathSeeds) {
                pathSeedSize += pathSeed.sourceLinks.size();

                // Select the pivot source link
                auto pivotSl = pathSeed.sourceLinks.at(0);

                // Find the ideal seed that contains the pivot
                auto idealSeedIt = std::find_if(
                    idealSeeds.begin(), 
                    idealSeeds.end(), 
                    [&](const auto& idealSeed) {
                        // Get the source links of the ideal seed
                        // and sort them by the geometry id
                        auto sls = idealSeed.sourceLinks;
                        std::sort(sls.begin(), sls.end(),
                            [&](const Acts::SourceLink& a, const Acts::SourceLink& b) {
                                auto sslA = a.get<SimpleSourceLink>();
                                auto sslB = b.get<SimpleSourceLink>();
                                std::size_t idA = static_cast<int>(sslA.geometryId().sensitive()/10 - 1);
                                std::size_t idB = static_cast<int>(sslB.geometryId().sensitive()/10 - 1);
                                return gOpt.staveZ.at(idA) < gOpt.staveZ.at(idB);
                        });

                        // Get the pivot source link
                        // of the path seed
                        auto ssl = pivotSl.get<SimpleSourceLink>();
                        
                        // Get the pivot source link
                        // of the ideal seed
                        Acts::SourceLink sl1 = sls.at(0);

                        auto ssl1 = sl1.get<SimpleSourceLink>();

                        // Compare the pivot source links
                        return ssl == ssl1;
                });
                // Pivot is not in the ideal seed
                // so the path seed is fake
                if (idealSeedIt == idealSeeds.end()) {
                    fakePivots++;
                    continue;
                }
                // Pivot is in the ideal seed
                matchedPivots++;
                
                // Compare the source links
                auto idealSeed = *idealSeedIt;
                bool seedMatched = true;
                for (int i = 0; i < idealSeed.sourceLinks.size(); i++) {
                    auto sl = idealSeed.sourceLinks.at(i);
                    auto ssl = sl.get<SimpleSourceLink>();
                    auto psIt = std::find_if(
                        pathSeed.sourceLinks.begin(), 
                        pathSeed.sourceLinks.end(), 
                        [&ssl](const Acts::SourceLink& psSl) {
                            auto psSsl = psSl.get<SimpleSourceLink>();
                            return ssl == psSsl;
                    });
                    if (psIt == pathSeed.sourceLinks.end()) {
                        lostSeeds++;
                        seedMatched = false;
                        break;
                    }
                }
                if (seedMatched) {
                    matchedSeeds++;
                }
            }

            std::cout << "Matched pivots: " << matchedPivots << " / " << idealPivots << std::endl;
            std::cout << "Fake pivots: " << fakePivots << std::endl;
            std::cout << "Lost pivots: " << lostPivots << std::endl;
            std::cout << "Matched seeds: " << matchedSeeds << " / " << idealPivots << std::endl;
            std::cout << "Lost seeds: " << lostSeeds << std::endl;
            std::cout << "Path seed size: " << pathSeedSize/pathSeeds.size() << std::endl;

            return ProcessCode::SUCCESS;
        }

    private:
        Config m_cfg;

        ReadDataHandle<Seeds> m_inputIdealSeeds
            {this, "InputIdealSeeds"};

        ReadDataHandle<Seeds> m_inputPathSeeds
            {this, "InputPathSeeds"};
};

int main() {
    // Set the log level
    Acts::Logging::Level logLevel = Acts::Logging::INFO;

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
    seqCfg.events = 1;
    seqCfg.numThreads = 1;
    seqCfg.trackFpes = false;
    Sequencer sequencer(seqCfg);

    // Add the sim data reader
    E320ROOTReader::E320ROOTSimDataReader::Config readerCfg = 
        E320ROOTReader::defaultSimConfig();
    readerCfg.dataCollection = "Measurements";
    std::string pathToDir = 
        "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/E320Pipeline_dataInRootFormat/Signal_E320lp_10.0";

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

    readerCfg.energyCuts = {0_GeV, 100_GeV};

    // Vertex position extent in the already rotated frame
    Acts::Extent vertexExtent;
    vertexExtent.set(
        Acts::BinningValue::binX, -100_mm, 100_mm);
    vertexExtent.set(
        Acts::BinningValue::binZ, -100_mm, 100_mm);
    vertexExtent.set(
        Acts::BinningValue::binY, -100_mm, 100_mm);

    readerCfg.vertexPosExtent = vertexExtent;

    // Add the reader to the sequencer
    sequencer.addReader(
        std::make_shared<E320ROOTReader::E320ROOTSimDataReader>(readerCfg, logLevel));

    // --------------------------------------------------------------
    // The path seeding setup
    SimpleSourceLink::SurfaceAccessor surfaceAccessor{*detector};

    auto pathSeederCfg = Acts::Experimental::PathSeeder<Grid>::Config();

    // Estimator of the IP and first hit
    // parameters of the track
    CsvLookupTableProvider::Config trackLookupCfg;

    trackLookupCfg.filePath = 
        "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_build/lookupTable.csv";

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
        geoId.setSensitive(i+1);
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
                {2, {300_um, 300_um}},
                {3, {350_um, 350_um}},
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
        .bins = std::make_pair(1000, 100),
    };
    gridConstructorCfg.surfaceAccessor.connect<
        &SimpleSourceLink::SurfaceAccessor::operator()>(
        &surfaceAccessor);

    auto gridConstructor = std::make_shared<E320TrackFinding::E320SourceLinkGridConstructor>(gridConstructorCfg); 

    // Create the path seeder algorithm
    auto seedingAlgoCfg = PathSeedingAlgorithm<Grid>::Config();
    seedingAlgoCfg.seeder = std::make_shared<Acts::Experimental::PathSeeder<Grid>>(pathSeederCfg);
    seedingAlgoCfg.sourceLinkGridConstructor = gridConstructor;
    seedingAlgoCfg.inputCollection = "Measurements";
    seedingAlgoCfg.outputCollection = "PathSeeds";

    sequencer.addAlgorithm(
        std::make_shared<PathSeedingAlgorithm<Grid>>(seedingAlgoCfg, logLevel));

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

    propOptions.maxSteps = 1000;

    auto options = Acts::KalmanFitterOptions(gctx, mctx, cctx, extensions,
        propOptions);

    // Reference surface for sampling the track at the IP
    double halfX = 10000;
    double halfY = 10000;

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
    TrackFitter<
        Propagator, 
        Trajectory, 
        TrackContainer>::Config fitterCfg{
            .inputCollection = "PathSeeds",
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
    // Ideal seeder for the truth comparison

    // Add the ideal seeder to the sequencer
    IdealSeeder::Config fullMatchingSeederCfg;

    fullMatchingSeederCfg.firstLayerIds = pathSeederCfg.firstLayerIds;
    fullMatchingSeederCfg.sourceLinkCalibrator.connect<
        &SimpleSourceLinkCoordinateCalibrator::operator()>(
        &sourceLinkCalibrator);
    fullMatchingSeederCfg.trackEstimator.connect<
        &CsvLookupTableProvider::operator()>(
        &trackLookup);

    auto fullMatchingSeeder = std::make_shared<IdealSeeder>(fullMatchingSeederCfg);

    IdealSeedingAlgorithm::Config idealSeederCfg{
        .seeder = fullMatchingSeeder,
        .inputCollection = "Measurements",
        .outputCollection = "IdealSeeds"};

    sequencer.addAlgorithm(
        std::make_shared<IdealSeedingAlgorithm>(idealSeederCfg, logLevel));

    // --------------------------------------------------------------
    // Event write out

    auto trackWriterCfg = ROOTFittedTrackWriter::Config{
        "SimTracks",
        "IdealSeeds",
        "fitted-tracks",
        "fitted-tracks.root",
        3,
        10
    };

    sequencer.addWriter(
        std::make_shared<ROOTFittedTrackWriter>(trackWriterCfg, logLevel));

    // // --------------------------------------------------------------
    // // Seed comparison

    // SeedComparer::Config comparerCfg{
        // .inputIdealCollection = "IdealSeeds",
        // .inputPathCollection = "PathSeeds"};

    // comparerCfg.surfaceAccessor.connect<
        // &SimpleSourceLink::SurfaceAccessor::operator()>(
        // &surfaceAccessor);

    // sequencer.addAlgorithm(
        // std::make_shared<SeedComparer>(comparerCfg, logLevel));

    return sequencer.run();
}
