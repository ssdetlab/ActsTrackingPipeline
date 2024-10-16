#pragma once

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"

#include "ActsLUXEPipeline/E320SourceLinkGrid.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Seeding/PathSeeder.hpp"

using namespace Acts::UnitLiterals;

class PathSeedingAlgorithm : public IAlgorithm {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// Path seeder
            std::shared_ptr<Acts::Experimental::PathSeeder> seeder;
            /// SourceLink grid
            std::shared_ptr<E320TrackFinding::E320SourceLinkGridConstructor> sourceLinkGridConstructor;
            /// Input source links
            std::string inputSourceLinks = "SourceLink";
            /// Output seeds
            std::string outputSeeds = "Seed";
        };

        /// @brief Constructor
        PathSeedingAlgorithm(Config config, Acts::Logging::Level level)
            : IAlgorithm("PathSeedingAlgorithm", level),
            m_cfg(std::move(config)) {
                m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
                m_outputSeeds.initialize(m_cfg.outputSeeds);
        }
        ~PathSeedingAlgorithm() = default;

        /// @brief The execute method        
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            // Get the input measurements
            // from the context
            auto input = m_inputSourceLinks(ctx);

            if (input.empty()) {
                m_outputSeeds(ctx, Seeds());
                return ProcessCode::SUCCESS;
            }

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

            auto gridLookup = m_cfg.sourceLinkGridConstructor->constructGrid(ctx.geoContext, input);

            std::vector<Acts::Experimental::PathSeeder::Seed> pathSeeds;
            m_cfg.seeder->getSeeds(ctx.geoContext, gridLookup, pathSeeds);

            Seeds outSeeds;
            auto me = 0.511 * Acts::UnitConstants::MeV;
            for (std::int32_t i = 0; i < pathSeeds.size(); i++) {
                const auto& seed = pathSeeds.at(i);
                Acts::Vector4 mPos4 = {seed.ipVertex.x(), seed.ipVertex.y(), seed.ipVertex.z(), 0};

                Acts::ActsScalar p = seed.ipP; 
                Acts::ActsScalar theta = std::acos(seed.ipDir.z()/p);
                Acts::ActsScalar phi = std::atan2(seed.ipDir.y(), seed.ipDir.x());

                Acts::CurvilinearTrackParameters ipParameters(
                    mPos4, phi, theta,
                    1_e / p, ipCov, 
                    Acts::ParticleHypothesis::electron());

                outSeeds.push_back(Seed{
                    seed.sourceLinks,
                    ipParameters,
                    i});
            }

            m_outputSeeds(ctx, std::move(outSeeds));
            
            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        Config m_cfg;

        ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks
            {this, "InputSourceLinks"};

        WriteDataHandle<Seeds> m_outputSeeds
            {this, "OutputSeeds"};
};
