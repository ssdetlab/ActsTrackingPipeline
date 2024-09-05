#pragma once

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"
#include "ActsLUXEPipeline/IdealSeeder.hpp"

#include "Acts/EventData/SourceLink.hpp"

/// @brief The ideal seeder for the LUXE simulation
/// takes the the SimMeasurements and converts them 
/// into seeds
class IdealSeedingAlgorithm : public IAlgorithm {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// Ideal seeder
            std::shared_ptr<IdealSeeder> seeder;
            /// The input collection
            std::string inputCollection = "SourceLink";
            /// The output collection
            std::string outputCollection = "Seed";
        };

        /// @brief Constructor
        IdealSeedingAlgorithm(Config config, Acts::Logging::Level level)
            : IAlgorithm("IdealSeedingAlgorithm", level),
            m_cfg(std::move(config)) {
                m_inputMeasurements.initialize(m_cfg.inputCollection);
                m_outputSeeds.initialize(m_cfg.outputCollection);
        }
        ~IdealSeedingAlgorithm() = default;

        /// @brief The execute method        
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            using namespace Acts::UnitLiterals;

            // Get the input measurements
            // from the context
            auto input = m_inputMeasurements(ctx);

            if (input.empty()) {
                m_outputSeeds(ctx, Seeds());
                return ProcessCode::SUCCESS;
            }

            // Sort the input by track id
            std::sort(input.begin(), input.end(),
                [](const auto& a, const auto& b) {
                    return a.trackId < b.trackId;
                }
            );

            Seeds seeds = m_cfg.seeder->getSeeds(ctx.geoContext, input);
            
            // std::cout << "IDEAL SEEDS SIZE: " << seeds.size() << std::endl;

            m_outputSeeds(ctx, std::move(seeds));

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;

        ReadDataHandle<SimMeasurements> m_inputMeasurements
            {this, "InputMeasurements"};

        WriteDataHandle<Seeds> m_outputSeeds
            {this, "OutputSeeds"};
};
