#pragma once

#include "Acts/EventData/SourceLink.hpp"

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"
#include "ActsLUXEPipeline/LUXEDataContainers.hpp"
#include "ActsLUXEPipeline/LUXEGeometryConstraints.hpp"

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
            /// Geometry options
            const LUXEGeometry::GeometryOptions& gOpt;
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

            // Sort the input by track id
            std::sort(input.begin(), input.end(),
                [](const auto& a, const auto& b) {
                    return a.trackId < b.trackId;
                }
            );

            // Create IP covariance matrix from 
            // reasonable standard deviations
            Acts::BoundVector ipStdDev;
            ipStdDev[Acts::eBoundLoc0] = 100_um;
            ipStdDev[Acts::eBoundLoc1] = 100_um;
            ipStdDev[Acts::eBoundTime] = 25_ns;
            ipStdDev[Acts::eBoundPhi] = 2_degree;
            ipStdDev[Acts::eBoundTheta] = 2_degree;
            ipStdDev[Acts::eBoundQOverP] = 1 / 100_GeV;
            Acts::BoundSquareMatrix ipCov = 
                ipStdDev.cwiseProduct(ipStdDev).asDiagonal();

            // Create the seeds
            LUXEDataContainer::Seeds seeds;
            std::vector<Acts::SourceLink> sourceLinks;
            for (auto it = input.begin(); it != input.end() - 1; ++it) {
                if (it->trackId == (it + 1)->trackId) {
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
                            continue;
                    }

                    // Ip parameter is the same for all hits
                    // with the same track id
                    Acts::CurvilinearTrackParameters ipParameters(
                        it->trueVertex, 
                        it->truthParameters[Acts::eBoundPhi],
                        it->truthParameters[Acts::eBoundTheta], 
                        -it->truthParameters[Acts::eBoundQOverP], 
                        ipCov,
                        Acts::ParticleHypothesis::electron());

                    // Add the seed to the list
                    // and reset the source links
                    seeds.push_back(LUXEDataContainer::Seed
                        {sourceLinks, ipParameters});
                    sourceLinks.clear();
                }
            }
            
            // Ip parameter is the same for all hits
            // with the same track id
            Acts::CurvilinearTrackParameters ipParameters(
                input.back().trueVertex, 
                input.back().truthParameters[Acts::eBoundPhi],
                input.back().truthParameters[Acts::eBoundTheta], 
                -input.back().truthParameters[Acts::eBoundQOverP], 
                ipCov,
                Acts::ParticleHypothesis::electron());

            // Add the last seed
            seeds.push_back(LUXEDataContainer::Seed
                {sourceLinks, ipParameters});

            m_outputSeeds(ctx, std::move(seeds));

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;

        ReadDataHandle<LUXEDataContainer::SimMeasurements> m_inputMeasurements
            {this, "InputMeasurements"};

        WriteDataHandle<LUXEDataContainer::Seeds> m_outputSeeds
            {this, "OutputSeeds"};
};
