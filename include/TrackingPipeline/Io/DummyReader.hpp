#pragma once

#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "Acts/EventData/SourceLink.hpp"

class DummyReader : public IReader {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// Number of dummy events to generate 
            size_t nEvents = 10;
            /// Collection with the measurement data
            std::string outputSourceLinks;
            /// Collection with the sim clusters data
            std::string outputSimClusters;
        };

        DummyReader(const DummyReader &) = delete;
        DummyReader(const DummyReader &&) = delete;

        /// Constructor
        /// @param config The Configuration struct
        /// @param level The log level
        DummyReader(const Config &config, Acts::Logging::Level level)
            : IReader(),
            m_cfg(config),
            m_logger(Acts::getDefaultLogger(name(), level)) {
                m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
                m_outputSimClusters.initialize(m_cfg.outputSimClusters);
        }

        /// Reader name() method
        virtual std::string name() const { return "DummyReader"; }
    
        /// Return the available events range.
        std::pair<std::size_t, std::size_t> 
            availableEvents() const override {
                return {0, m_cfg.nEvents};
        }

        /// Read out data from the input stream
        ProcessCode read(const AlgorithmContext &context) override {
            // Create empty measurements
            std::vector<Acts::SourceLink> sourceLinks;
            SimClusters clusters;

            m_outputSourceLinks(context, std::move(sourceLinks));
            m_outputSimClusters(context, std::move(clusters));

            // Return success flag
            return ProcessCode::SUCCESS;
        }

    private:
        /// Private access to the logging instance
        const Acts::Logger &logger() const { return *m_logger; }

        /// The config class
        Config m_cfg;

        /// WriteDataHandle for the source links
        WriteDataHandle<std::vector<Acts::SourceLink>> m_outputSourceLinks
            {this, "OutputObsData"};
        
        /// WriteDataHandle for the sim clusters
        WriteDataHandle<SimClusters> m_outputSimClusters
            {this, "OutputSimClusters"};

        std::unique_ptr<const Acts::Logger> m_logger;
};
