#pragma once

#include "Acts/EventData/SourceLink.hpp"

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"
#include "ActsLUXEPipeline/LUXEDataContainers.hpp"

#include <fstream>

/// @brief The ideal seeder for the LUXE simulation
/// takes the the SimMeasurements and converts them 
/// into seeds
class IdealSeeder : public IAlgorithm {
    public:
        /// @brief The nested configuration struct
        struct Config {
            std::string inputCollection = "SourceLink";
            std::string outputCollection = "Seed";
            std::string lookupTableDir = "LUXEPipeline/build";
            std::uint32_t minHits = 4;
            std::uint32_t maxHits = 8;
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
            // Get the input measurements
            // from the context
            auto input = m_inputMeasurements(ctx);

            // Sort the input by track id
            std::sort(input.begin(), input.end(),
                [](const auto& a, const auto& b) {
                    return a.trackId < b.trackId;
                }
            );

            // Create the seeds
            LUXEDataContainer::Seeds seeds;
            std::vector<Acts::SourceLink> sourceLinks{input.at(0).sourceLink};
            // Ip parameter is the same for all hits
            // with the same track id

            std::unordered_map <float, float> EX1_map =
                    readLookup(m_cfg.lookupTableDir+"EX1_lookup_table.txt");
            std::unordered_map <float, float> X1X4_map =
                    readLookup(m_cfg.lookupTableDir+"X1X4_lookup_table.txt");
            std::unordered_map <float, float> Z1Z4_map =
                    readLookup(m_cfg.lookupTableDir+"Z1Z4_lookup_table.txt");

            Acts::Vector3 ipParameters{
                input.at(0).truthParameters[Acts::eBoundQOverP],
                input.at(0).truthParameters[Acts::eBoundPhi],
                input.at(0).truthParameters[Acts::eBoundTheta]};
            for (auto it = input.begin() + 1; it != input.end(); ++it) {
                if (it->trackId == (it - 1)->trackId) {
                    // Add the source link to the list
                    // if the hit is from the same track
                    sourceLinks.push_back(it->sourceLink);
                }
                else {
                    if (sourceLinks.size() < m_cfg.minHits || 
                        sourceLinks.size() > m_cfg.maxHits) {
                            // If the seed does not have the right size
                            // skip it
                            sourceLinks.clear();
                            sourceLinks.push_back(it->sourceLink);
                            ipParameters.x() = it->truthParameters[Acts::eBoundQOverP],
                            ipParameters.y() = it->truthParameters[Acts::eBoundPhi],
                            ipParameters.z() = it->truthParameters[Acts::eBoundTheta];
                            continue;
                    }
                    // Add the seed to the list
                    // and reset the source links
                    seeds.push_back(LUXEDataContainer::Seed
                        {sourceLinks, ipParameters});
                    sourceLinks.clear();
                    sourceLinks.push_back(it->sourceLink);
                    ipParameters.x() = it->truthParameters[Acts::eBoundQOverP],
                    ipParameters.y() = it->truthParameters[Acts::eBoundPhi],
                    ipParameters.z() = it->truthParameters[Acts::eBoundTheta];
                }
            }
            
            // Add the last seed
            seeds.push_back(LUXEDataContainer::Seed
                {sourceLinks, ipParameters});

            m_outputSeeds(ctx, std::move(seeds));

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

        const std::unordered_map <float, float> readLookup(std::string file) const {
            std::ifstream lookupFile(&file[0]);
            std::unordered_map <float, float> lookupTable;
            float x, y;
            while (lookupFile >> x >> y) {
                lookupTable[x] = y;
            }
            lookupFile.close();
            return lookupTable;
        }
    private:
        Config m_cfg;

        ReadDataHandle<LUXEDataContainer::SimMeasurements> m_inputMeasurements
            {this, "InputMeasurements"};

        WriteDataHandle<LUXEDataContainer::Seeds> m_outputSeeds
            {this, "OutputSeeds"};
};
