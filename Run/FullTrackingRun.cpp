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
#include "ActsLUXEPipeline/IdealSeedingAlgorithm.hpp"
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

            std::cout << "Ideal seeds size: " << idealSeeds.size() << std::endl;
            std::cout << "Path seeds size: " << pathSeeds.size() << std::endl;

            if (idealSeeds.empty() || pathSeeds.empty()) {
                return ProcessCode::SUCCESS;
            }

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
            int lostSeeds23 = 0;
            int lostSeeds34 = 0;

            double seeds23 = 0;
            double seeds34 = 0;

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
                if (idealSeed.ipParameters.absoluteMomentum() > 2 && 
                    idealSeed.ipParameters.absoluteMomentum() < 3) {
                        seeds23++;
                }
                else if (idealSeed.ipParameters.absoluteMomentum() > 3 &&
                    idealSeed.ipParameters.absoluteMomentum() < 4) {
                        seeds34++;
                }
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
                        double E = idealSeed.ipParameters.absoluteMomentum();
                        if (E > 2 && E < 3) {
                            lostSeeds23++;
                        }
                        else if (E > 3 && E < 4) {
                            lostSeeds34++;
                        }
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
            std::cout << "Lost seeds 23: " << 1 - lostSeeds23 / seeds23 << std::endl;
            std::cout << "Lost seeds 34: " << 1 - lostSeeds34 / seeds34 << std::endl;
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
    seqCfg.events = 100;
    seqCfg.numThreads = -1;
    seqCfg.trackFpes = false;
    Sequencer sequencer(seqCfg);

    // Add the sim data reader
    E320ROOTReader::E320ROOTSimSplitDataReader::Config readerCfg = 
        E320ROOTReader::defaultSimSplitConfig();
    readerCfg.dataCollection = "SimMeasurements";
    std::string pathToDir = 
        "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/E320Pipeline_dataInRootFormat/Signal_E320lp_10.0_split";

    // Get the paths to the files in the directory
    for (const auto & entry : std::filesystem::directory_iterator(pathToDir)) {
        std::string pathToFile = entry.path();
        readerCfg.filePaths.push_back(pathToFile);
        std::cout << pathToFile << std::endl;
    }

    std::cout << "Number of files: " << readerCfg.filePaths.size() << std::endl;

    // // The events are not sorted in the directory
    // // but we need to process them in order
    // std::sort(readerCfg.filePaths.begin(), readerCfg.filePaths.end(),
        // [] (const std::string& a, const std::string& b) {
            // std::size_t idxRootA = a.find_last_of('.');
            // std::size_t idxEventA = a.find_last_of('t', idxRootA);
            // std::string eventSubstrA = a.substr(idxEventA + 1, idxRootA - idxEventA);
            
            // std::size_t idxRootB = b.find_last_of('.');
            // std::size_t idxEventB = b.find_last_of('t', idxRootB);
            // std::string eventSubstrB = b.substr(idxEventB + 1, idxRootB - idxEventB);

            // return std::stoul(eventSubstrA) < std::stoul(eventSubstrB);
        // }
    // );

    readerCfg.energyCuts = {2_GeV, 4_GeV};

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
        std::make_shared<E320ROOTReader::E320ROOTSimSplitDataReader>(readerCfg, logLevel));

    // --------------------------------------------------------------
    // Backgorund generation
    auto noiseGenerator = std::make_shared<UniformNoiseGenerator>();
    noiseGenerator->numberOfHits = 500;

    NoiseEmbeddingAlgorithm::Config noiseCfg {
        .noiseGenerator = noiseGenerator,
        .detector = detector.get(),
        .inputCollection = "SimMeasurements",
        .outputCollection = "Measurements"
    };

    sequencer.addAlgorithm(
        std::make_shared<NoiseEmbeddingAlgorithm>(noiseCfg, logLevel));

    // --------------------------------------------------------------
    // The path seeding setup
    SimpleSourceLink::SurfaceAccessor surfaceAccessor{*detector};

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

    pathSeederCfg.minE = 1.8_GeV;
    pathSeederCfg.maxE = 4.2_GeV;

    // Path width provider
    std::map<std::int32_t, std::pair<
        Acts::ActsScalar,Acts::ActsScalar>> 
            pathWidths = {
                {0, {100_um, 100_um}},
                {1, {200_um, 200_um}},
                {2, {450_um, 450_um}},
                {3, {700_um, 700_um}},
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
        .bins = std::make_pair(1000, 10),
    };
    gridConstructorCfg.surfaceAccessor.connect<
        &SimpleSourceLink::SurfaceAccessor::operator()>(
        &surfaceAccessor);

    auto gridConstructor = std::make_shared<E320TrackFinding::E320SourceLinkGridConstructor>(gridConstructorCfg); 

    // Create the path seeder algorithm
    auto seedingAlgoCfg = PathSeedingAlgorithm::Config();
    seedingAlgoCfg.seeder = std::make_shared<Acts::Experimental::PathSeeder>(pathSeederCfg);
    seedingAlgoCfg.sourceLinkGridConstructor = gridConstructor;
    seedingAlgoCfg.inputCollection = "Measurements";
    seedingAlgoCfg.outputCollection = "PathSeeds";

    sequencer.addAlgorithm(
        std::make_shared<PathSeedingAlgorithm>(seedingAlgoCfg, logLevel));

    // --------------------------------------------------------------
    // Track finding
    // TryAllTrackFindingAlgorithm::Config trackFindingCfg;
    // trackFindingCfg.inputCollection = "PathSeeds";
    // trackFindingCfg.outputCollection = "TrackCandidates";
    // trackFindingCfg.minSourceLinks = 3;

    // sequencer.addAlgorithm(
        // std::make_shared<TryAllTrackFindingAlgorithm>(trackFindingCfg, logLevel));

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
                            {10}, 
                            {10u}
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
    trackFindingCfg.inputCollection = "PathSeeds";
    trackFindingCfg.outputCollection = "TrackCandidates";
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
    // Ideal seeder for the truth comparison

    // Add the ideal seeder to the sequencer
    IdealSeeder::Config fullMatchingSeederCfg;

    fullMatchingSeederCfg.minSourceLinks = 4;
    fullMatchingSeederCfg.maxSourceLinks = 4;

    auto firstTrackingVolume = detector->findDetectorVolume("layer0");
    for (const auto& s : firstTrackingVolume->surfaces()) {
        fullMatchingSeederCfg.firstLayerIds.push_back(s->geometryId());
    }
    // fullMatchingSeederCfg.sourceLinkCalibrator.connect<
        // &SimpleSourceLinkCoordinateCalibrator::operator()>(
        // &sourceLinkCalibrator);
    // fullMatchingSeederCfg.trackEstimator.connect<
        // &CsvLookupTableProvider::operator()>(
        // &trackLookup);

    auto fullMatchingSeeder = std::make_shared<IdealSeeder>(fullMatchingSeederCfg);

    IdealSeedingAlgorithm::Config idealSeederCfg{
        .seeder = fullMatchingSeeder,
        .inputCollection = "Measurements",
        .outputCollection = "IdealSeeds"};

    sequencer.addAlgorithm(
        std::make_shared<IdealSeedingAlgorithm>(idealSeederCfg, logLevel));
        
    // --------------------------------------------------------------
    // Event write out

    auto trackWriterCfg = ROOTFittedTrackWriter::Config();
    trackWriterCfg.surfaceAccessor.connect<
        &SimpleSourceLink::SurfaceAccessor::operator()>(
            &surfaceAccessor);

    trackWriterCfg.inputTrackCollection = "Tracks";
    trackWriterCfg.inputSeedCollection = "IdealSeeds";
    trackWriterCfg.treeName = "fitted-tracks";
    trackWriterCfg.filePath = "fitted-tracks-signal-bkg.root";

    sequencer.addWriter(
        std::make_shared<ROOTFittedTrackWriter>(trackWriterCfg, logLevel));

    // // --------------------------------------------------------------
    // // Phoenix write out
    // PhoenixTrackWriter::Config phoenixWriterCfg;
    // phoenixWriterCfg.inputTrackCollection = "Tracks";
    // phoenixWriterCfg.fileName = "test-tracks";

    // sequencer.addWriter(
        // std::make_shared<PhoenixTrackWriter>(phoenixWriterCfg, logLevel));

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
