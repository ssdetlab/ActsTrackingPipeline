#include "ActsLUXEPipeline/LUXEROOTDataReader.hpp"
#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXEIdealSeeder.hpp"
#include "ActsLUXEPipeline/Sequencer.hpp"

#include <filesystem>

class DummyAlgo : public IAlgorithm {
    public:
        struct Config {
            std::string inputCollection = "SimSeeds";
            const Acts::Experimental::Detector& detector;
        };

        DummyAlgo(Config config, Acts::Logging::Level level)
            : IAlgorithm("Dummy", level),
            m_cfg(std::move(config)) {
                m_inputCollection.initialize(m_cfg.inputCollection);
        }
        ~DummyAlgo() = default;

        ProcessCode execute(const AlgorithmContext& ctx) const override {
            auto input = m_inputCollection(ctx);

            auto surfaceAccessor = 
                SimpleSourceLink::SurfaceAccessor{m_cfg.detector};

            for (const auto& seed : input) {
                ACTS_INFO("Seed with " << seed.sourceLinks.size() << " source links");
                for (const auto& sl : seed.sourceLinks) {
                    auto surface = surfaceAccessor(sl);
                    assert(surface != nullptr);
                    ACTS_INFO("Surface: " << surface->geometryId());
                }
            }

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;

        ReadDataHandle<LUXEDataContainer::Seeds> m_inputCollection{
            this, "InputSeeds"};

};

int main() {
    // Set the log level
    Acts::Logging::Level logLevel = Acts::Logging::VERBOSE;

    // Dummy context and options
    Acts::GeometryContext gctx;
    LUXEGeometry::GeometryOptions gOpt;

    // Set the path to the gdml file
    // and the names of the volumes to be converted
    std::string gdmlPath = 
        "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_gdmls/lxgeomdump_stave_positron.gdml";
    std::vector<std::string> names{"OPPPSensitive"};

    // Build the LUXE detector
    auto positronArmBpr = 
        LUXEGeometry::makeBlueprintPositron(gdmlPath, names, gOpt);
    auto detector =
        LUXEGeometry::buildLUXEDetector(std::move(positronArmBpr), gctx, gOpt);

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
    IdealSeeder::Config seederCfg;
    seederCfg.inputCollection = "SimMeasurements";
    seederCfg.outputCollection = "SimSeeds";
    sequencer.addAlgorithm(
        std::make_shared<IdealSeeder>(seederCfg, logLevel));

    // Add the dummy algorithm to the sequencer
    DummyAlgo::Config dummyCfg{
        .inputCollection = "SimSeeds",
        .detector = *detector};
    
    dummyCfg.inputCollection = "SimSeeds";

    sequencer.addAlgorithm(
        std::make_shared<DummyAlgo>(dummyCfg, logLevel));

    // Run all configured algorithms and return the appropriate status.
    return sequencer.run();
}