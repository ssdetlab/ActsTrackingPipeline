#pragma once

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"

#include "Acts/EventData/SourceLink.hpp"

/// @brief The ideal seeder for the LUXE simulation
/// takes the the SimMeasurements and converts them 
/// into seeds
class IdealSeeder : public IAlgorithm {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// The input collection
            std::string inputCollection = "SourceLink";
            /// The output collection
            std::string outputCollection = "Seed";
            /// The minimum number of hits
            /// for a seed to be created
            std::uint32_t minHits = 4;
            /// The maximum number of hits
            /// for a seed to be created
            std::uint32_t maxHits = 4;
        };

        /// @brief Constructor
        IdealSeeder(Config config, Acts::Logging::Level level)
            : IAlgorithm("IdealSeeder", level),
            m_cfg(std::move(config)) {
                m_inputMeasurements.initialize(m_cfg.inputCollection);
                m_outputSeeds.initialize(m_cfg.outputCollection);
        }
        ~IdealSeeder() = default;

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

            // Create the seeds
            Seeds seeds;
            std::vector<Acts::SourceLink> sourceLinks;

            // Insert the first source link
            sourceLinks.push_back(input.front().sourceLink);
            for (auto it = input.begin() + 1; it != input.end(); ++it) {
                if (it->trackId == (it - 1)->trackId && (it + 1) != input.end()) {
                    // Add the source link to the list
                    // if the hit is from the same track
                    sourceLinks.push_back(it->sourceLink);
                }
                else {
                    if ((it + 1) == input.end()) {
                        // Add the last source link
                        sourceLinks.push_back(it->sourceLink);
                    }
                    if (sourceLinks.size() < m_cfg.minHits || 
                        sourceLinks.size() > m_cfg.maxHits) {
                            // If the seed does not have the right size
                            // skip it
                            sourceLinks.clear();
                            continue;
                    }

                    // Ip parameter is the same for all hits
                    // with the same track id
                    Acts::CurvilinearTrackParameters ipParameters =
                        (it - 1)->ipParameters;

                    // Add the seed to the list
                    // and reset the source links
                    seeds.push_back(Seed
                        {sourceLinks, ipParameters, (it - 1)->trackId});
                    sourceLinks.clear();
                    sourceLinks.push_back(it->sourceLink);
                }
            }
            
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
