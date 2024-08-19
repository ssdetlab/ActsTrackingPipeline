#include "ActsLUXEPipeline/Sequencer.hpp"
#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/BinnedMagneticField.hpp"
#include "ActsLUXEPipeline/ConstantBoundedField.hpp"
#include "ActsLUXEPipeline/MeasurementsCreator.hpp"
#include "ActsLUXEPipeline/AlgorithmContext.hpp"
#include "ActsLUXEPipeline/ROOTLookupDataWriter.hpp"
#include "ActsLUXEPipeline/Generators.hpp"
#include "ActsLUXEPipeline/LUXESourceLinkBinner.hpp"
#include "ActsLUXEPipeline/PathSeeder.hpp"
#include "ActsLUXEPipeline/IdealSeeder.hpp"
#include "ActsLUXEPipeline/CsvLookupTableProvider.hpp"
#include "ActsLUXEPipeline/ForwardOrderedIntersectionFinder.hpp"
#include "ActsLUXEPipeline/LUXEPathWidthProvider.hpp"
#include "ActsLUXEPipeline/TrackFitter.hpp"
#include "ActsLUXEPipeline/ROOTFittedTrackWriter.hpp"

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"

#include <filesystem>

using namespace Acts::UnitLiterals;

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

using TrackParameters = Acts::CurvilinearTrackParameters;

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

            std::cout << "Ideal seeds = " << idealSeeds.size() << std::endl;
            std::cout << "Path seeds = " << pathSeeds.size() << std::endl;

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
            int lostSeeds = 0;

            // Normalization
            int idealPivots = idealSeeds.size();

            LUXEGeometry::GeometryOptions gOpt;

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
                if (id != 0 && id != 1) {
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

                        // Compare the pivot source links
                        return ssl == ssl1;
                });
                // Pivot is not in the ideal seed
                // so the path seed is fake
                if (pathSeedIt == pathSeeds.end()) {
                    lostPivots++;
                    continue;
                }
            }

            // Iterate over the path seeds
            for (auto pathSeed : pathSeeds) {
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

            return ProcessCode::SUCCESS;
        }

    private:
        Config m_cfg;

        ReadDataHandle<Seeds> m_inputIdealSeeds
            {this, "InputIdealSeeds"};

        ReadDataHandle<Seeds> m_inputPathSeeds
            {this, "InputPathSeeds"};
};

/// @brief Run the propagation through 
/// a uniform energy spectrum and record the
/// measurements. Run two seeders and compare
/// the results
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

    std::string materialPath = "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_material/material_uniform_binned_16X_8Y.json";

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
        Acts::BinningValue::binX, 
        gOpt.dipoleTranslation[0] - gOpt.dipoleBounds[0] + gOpt.constantFieldDelta[0],
        gOpt.dipoleTranslation[0] + gOpt.dipoleBounds[0] - gOpt.constantFieldDelta[0]);
    dipoleExtent.set(
        Acts::BinningValue::binZ,
        gOpt.dipoleTranslation[1] - gOpt.dipoleBounds[1] + gOpt.constantFieldDelta[1],
        gOpt.dipoleTranslation[1] + gOpt.dipoleBounds[1] - gOpt.constantFieldDelta[1]);
    dipoleExtent.set(
        Acts::BinningValue::binY,
        gOpt.dipoleTranslation[2] - gOpt.dipoleBounds[2] + gOpt.constantFieldDelta[2],
        gOpt.dipoleTranslation[2] + gOpt.dipoleBounds[2] - gOpt.constantFieldDelta[2]);

    auto field = std::make_shared<ConstantBoundedField>(
        Acts::Vector3(0., 0., -1.2_T),
        dipoleExtent);

    // --------------------------------------------------------------
    // Event creation 

    // Setup the sequencer
    Sequencer::Config seqCfg;
    seqCfg.events = 10000;
    seqCfg.numThreads = 1;
    seqCfg.trackFpes = false;
    Sequencer sequencer(seqCfg);

    // Setup the measurements creator
    Acts::Experimental::DetectorNavigator::Config navCfg;
    navCfg.detector = detector.get();
    navCfg.resolvePassive = false;
    navCfg.resolveMaterial = true;
    navCfg.resolveSensitive = true;

    Acts::Experimental::DetectorNavigator navigator(
        navCfg, 
        Acts::getDefaultLogger(
            "DetectorNavigator", 
            logLevel));
    Acts::EigenStepper<> stepper(field);

    auto propagator = 
        Propagator(std::move(stepper), std::move(navigator));

    MeasurementsCreator::Config mcCfg;
    mcCfg.outputCollection = "Measurements";
    mcCfg.vertexGenerator = std::make_shared<StationaryVertexGenerator>();
    mcCfg.momentumGenerator = std::make_shared<LUXESimParticle::RangedUniformMomentumGenerator>();
    mcCfg.randomNumberSvc = std::make_shared<RandomNumbers>(RandomNumbers::Config());
    mcCfg.nTracks = 1;

    sequencer.addAlgorithm(
        std::make_shared<MeasurementsCreator>(
                propagator, mcCfg, logLevel));

    // --------------------------------------------------------------
    // The path seeding setup

    // Setup the source link binner
    SimpleSourceLink::SurfaceAccessor surfaceAccessor{*detector};

    LUXETrackFinding::LUXESourceLinkBinner::Config sourceLinkBinnerCfg{
        .gOpt = gOpt,
        .bins = std::make_pair(300, 10),
    };
    sourceLinkBinnerCfg.surfaceAccessor.connect<
        &SimpleSourceLink::SurfaceAccessor::operator()>(
        &surfaceAccessor);

    auto sourceLinkBinner = 
        std::make_shared<LUXETrackFinding::LUXESourceLinkBinner>(sourceLinkBinnerCfg);

    // Setup the path seeder
    PathSeeder::Config pathSeederCfg;
    pathSeederCfg.inputCollection = "Measurements";
    pathSeederCfg.outputCollection = "PathSeeds";
    pathSeederCfg.minHits = 3;
    pathSeederCfg.maxHits = 1e8;
    pathSeederCfg.sourceLinkBinner = sourceLinkBinner;
    pathSeederCfg.surfaceAccessor.connect<
        &SimpleSourceLink::SurfaceAccessor::operator()>(
        &surfaceAccessor);

    // Lookup table providers
    auto exfLookupCfg = CsvLookupTableProvider::Config{
        .filePath = "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_lookups/csvs/200k_ranged_uniform/XFirstE_lookup_table.csv"
    };
    auto xfxLookupCfg = CsvLookupTableProvider::Config{
        .filePath = "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_lookups/csvs/200k_ranged_uniform/XFirstXLast_lookup_table.csv"
    };
    auto xfyLookupCfg = CsvLookupTableProvider::Config{
        .filePath = "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_lookups/csvs/200k_ranged_uniform/XFirstYLast_lookup_table.csv"
    };
    auto zfLookupCfg = CsvLookupTableProvider::Config{
        .filePath = "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_lookups/csvs/200k_ranged_uniform/ZFirstZLast_lookup_table.csv"
    };

    pathSeederCfg.EXFirstLookupProvider =
        std::make_shared<CsvLookupTableProvider>(exfLookupCfg);
    pathSeederCfg.XFirstXLastLookupProvider =
        std::make_shared<CsvLookupTableProvider>(xfxLookupCfg);
    pathSeederCfg.XFirstYLastLookupProvider =
        std::make_shared<CsvLookupTableProvider>(xfyLookupCfg);
    pathSeederCfg.ZFirstZLastLookupProvider =
        std::make_shared<CsvLookupTableProvider>(zfLookupCfg);

    // Path width provider
    std::map<std::int32_t, std::pair<
        Acts::ActsScalar,Acts::ActsScalar>> 
            pathWidths = {
                {0, {250_um, 250_um}},
                {1, {250_um, 250_um}},
                {2, {300_um, 300_um}},
                {3, {300_um, 300_um}},
                {4, {300_um, 300_um}},
                {5, {400_um, 400_um}},
                {6, {400_um, 400_um}},
                {7, {400_um, 400_um}}
    };

    auto pathWidthProvider = std::make_shared<LUXETrackFinding::LUXEPathWidthProvider>(
        gOpt,
        pathWidths);

    pathSeederCfg.pathWidthProvider = pathWidthProvider;
    
    // Intersection finder
    auto intesectionFinder = std::make_shared<ForwardOrderedIntersectionFinder>();

    // Combine layers into surfaces for the intersection finder
    std::vector<std::shared_ptr<Acts::Surface>> surfacePtrs;
    for (int i = 0; i < gOpt.staveZ.size(); i++) {
        double halfX = 
            ((gOpt.chipXEven.at(8) + gOpt.chipSizeX/2) - 
            (gOpt.chipXEven.at(0) - gOpt.chipSizeX/2))/2;

        double halfY = gOpt.chipSizeY/2;

        double centerX = (i % 2 == 0) ?
            ((gOpt.chipXEven.at(8) + gOpt.chipSizeX/2) + 
            (gOpt.chipXEven.at(0) - gOpt.chipSizeX/2))/2 :
            ((gOpt.chipXOdd.at(8) + gOpt.chipSizeX/2) +
            (gOpt.chipXOdd.at(0) - gOpt.chipSizeX/2))/2;

        Acts::Transform3 transform(
            Acts::Translation3(Acts::Vector3(centerX, gOpt.staveZ.at(i), 0)) *
            gOpt.actsToWorldRotation.inverse());

        auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
            transform,
            std::make_shared<Acts::RectangleBounds>(
                halfX, halfY));

        Acts::GeometryIdentifier geoId;
        geoId.setSensitive(i+1);
        surface->assignGeometryId(std::move(geoId));
        surfacePtrs.push_back(surface);
    }
    std::vector<const Acts::Surface*> surfaces;
    for (auto surface : surfacePtrs) {
        surfaces.push_back(surface.get());
    }

    intesectionFinder->m_surfaces = std::move(surfaces);
    intesectionFinder->tol = 
        (gOpt.chipXOdd.at(1) - gOpt.chipSizeX/2) - (gOpt.chipXOdd.at(0) + gOpt.chipSizeX/2) + 1_mm;
    pathSeederCfg.intersectionFinder = intesectionFinder;

    // Extent in already rotated frame
    Acts::Extent firstLayerExtent;
    firstLayerExtent.set(
        Acts::BinningValue::binX, 
        gOpt.chipXEven.at(0) - gOpt.chipSizeX,
        gOpt.chipXOdd.at(8) + gOpt.chipSizeX);
    firstLayerExtent.set(
        Acts::BinningValue::binZ,
        -gOpt.chipY - gOpt.chipSizeY,
        -gOpt.chipY + gOpt.chipSizeY);
    firstLayerExtent.set(
        Acts::BinningValue::binY,
        gOpt.layerZPositions.at(0) - gOpt.layerBounds.at(2) - 1_mm,
        gOpt.layerZPositions.at(0) + gOpt.layerBounds.at(2) + 1_mm);

    pathSeederCfg.firstLayerExtent = firstLayerExtent;

    // Add the path seeder to the sequencer
    sequencer.addAlgorithm(
        std::make_shared<PathSeeder>(pathSeederCfg, logLevel));

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

    propOptions.maxSteps = 200;

    auto options = Acts::KalmanFitterOptions(gctx, mctx, cctx, extensions,
        propOptions);

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

    // // Add the ideal seeder to the sequencer
    // IdealSeeder::Config seederCfg{
        // .inputCollection = "Measurements",
        // .outputCollection = "IdealSeeds",
        // .minHits = 3,
        // .maxHits = 10,
    // };
    // sequencer.addAlgorithm(
        // std::make_shared<IdealSeeder>(seederCfg, logLevel));

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
