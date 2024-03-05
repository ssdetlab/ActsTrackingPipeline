#include "ActsLUXEPipeline/LUXEROOTDataReader.hpp"
#include "ActsLUXEPipeline/Sequencer.hpp"

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"

#include <filesystem>

class dummyAlgorithm : public IAlgorithm {
    public:
        struct Config {
            std::string inputSourceLinks = "SourceLink";
        };

        dummyAlgorithm(Config config, Acts::Logging::Level level)
            : IAlgorithm("dummyAlgorithm", level),
            m_cfg(std::move(config)) {
            m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
        }
        ~dummyAlgorithm() = default;

        ProcessCode execute(const AlgorithmContext& ctx) const override {
            const auto& sourceLinks = m_inputSourceLinks(ctx);
            
            for (const auto& sourceLink : sourceLinks.sourceLinks) {
                auto idx = sourceLink.get<SimpleSourceLink>();
                std::cout << "SourceLink: " << idx.geometryId() << " " << idx.index() << std::endl;
            }

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;

        ReadDataHandle<LUXEROOTReader::SimMeasurements> m_inputSourceLinks{
            this, "InputSourceLinks"};
};

/// @brief Run the propagation through 
/// a uniform energy spectrum and record the
/// energy vs position histograms for each layer
int main() {
    Acts::Logging::Level logLevel = Acts::Logging::VERBOSE;

    // setup the sequencer first w/ config derived from options
    Sequencer::Config seqCfg;
    seqCfg.events = 10;
    seqCfg.numThreads = -1;
    Sequencer sequencer(seqCfg);

    LUXEROOTReader::LUXEROOTSimDataReader::Config readerCfg 
        = LUXEROOTReader::defaultSimConfig();
    readerCfg.dataCollection = "SourceLink";
    std::string pathToDir = "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_dataInRootFormat/SignalNextTrial_e1gpc_10.0_1";

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

    readerCfg.filePaths = std::vector<std::string>(
        readerCfg.filePaths.begin(), readerCfg.filePaths.begin() + 72); 

    sequencer.addReader(
        std::make_shared<LUXEROOTReader::LUXEROOTSimDataReader>(readerCfg, logLevel));

    dummyAlgorithm::Config algCfg;
    algCfg.inputSourceLinks = "SourceLink";
    sequencer.addAlgorithm(
        std::make_shared<dummyAlgorithm>(algCfg, logLevel));

    // Run all configured algorithms and return the appropriate status.
    return sequencer.run();
}