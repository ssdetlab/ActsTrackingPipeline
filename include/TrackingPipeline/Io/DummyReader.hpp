#pragma once

#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

#include <cstddef>

using namespace Acts::UnitLiterals;

class DummyReader : public IReader {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// Output source links
            std::string outputSourceLinks;
            /// Output sim clusters
            std::string outputSimClusters;
            /// Number of events
            std::size_t nEvents;
        };

        /// @brief Constructor
        DummyReader(
            const Config& config)
            : IReader(),
            m_cfg(config) {
                m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
                m_outputSimClusters.initialize(m_cfg.outputSimClusters);
        }

        /// @brief The execute method        
        ProcessCode read(const AlgorithmContext& ctx) override {
            std::vector<Acts::SourceLink> sourceLinks{};
            SimClusters clusters{};

            m_outputSourceLinks(ctx, std::move(sourceLinks));
            m_outputSimClusters(ctx, std::move(clusters));

            return ProcessCode::SUCCESS;
        }

        /// @brief Provide range of available events or [0, SIZE_MAX) if undefined.
        std::pair<std::size_t, std::size_t> availableEvents() const override {
            return {0, m_cfg.nEvents};
        }

        /// @brief Reader name.
        std::string name() const override { return "DummyReader"; }

        /// @brief Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        Config m_cfg;

        WriteDataHandle<std::vector<Acts::SourceLink>> m_outputSourceLinks{
            this, "OutputSourceLinks"};

        WriteDataHandle<SimClusters> m_outputSimClusters{
            this, "OutputSimClusters"};
};
