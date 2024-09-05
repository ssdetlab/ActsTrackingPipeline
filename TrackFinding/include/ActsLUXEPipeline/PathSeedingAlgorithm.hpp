#pragma once

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"

#include "ActsLUXEPipeline/E320SourceLinkGrid.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Seeding/PathSeeder.hpp"

struct TryAllTracks {
    void tryAllTracks(
        std::shared_ptr<Acts::Experimental::SeedTree::Node> root,
        std::vector<Acts::SourceLink> track) {
            auto parentPars = root->m_sourceLink.get<SimpleSourceLink>().parameters;

            if (root->children.size() == 0) {
                // std::cout << "Parent: (" << parentPars[0] << " " << parentPars[1]
                // << ") -----> NO CHILDREN" << std::endl;

                track.push_back(root->m_sourceLink);
                tracks.push_back(track);
            }

            track.push_back(root->m_sourceLink);
            for (auto& child : root->children) {
                auto childPars = child->m_sourceLink.get<SimpleSourceLink>().parameters;
                // std::cout << "Parent: (" << parentPars[0] << " " << parentPars[1]
                // << ") -----> THE CHILD: (" << childPars[0] << " " << childPars[1] << ")" << std::endl;

                tryAllTracks(child, track);
            }    
    }

    std::vector<std::vector<Acts::SourceLink>> tracks;
};

template <typename grid_t>
class PathSeedingAlgorithm : public IAlgorithm {
    public:
        using GridType = grid_t;

        /// @brief The nested configuration struct
        struct Config {
            /// Path seeder
            std::shared_ptr<Acts::Experimental::PathSeeder<GridType>> seeder;
            /// SourceLink grid
            std::shared_ptr<E320TrackFinding::E320SourceLinkGridConstructor> sourceLinkGridConstructor;
            /// The input collection
            std::string inputCollection = "SourceLink";
            /// The output collection
            std::string outputCollection = "Seed";
        };

        /// @brief Constructor
        PathSeedingAlgorithm(Config config, Acts::Logging::Level level)
            : IAlgorithm("PathSeedingAlgorithm", level),
            m_cfg(std::move(config)) {
                m_inputMeasurements.initialize(m_cfg.inputCollection);
                m_outputSeeds.initialize(m_cfg.outputCollection);
        }
        ~PathSeedingAlgorithm() = default;

        /// @brief The execute method        
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            // Get the input measurements
            // from the context
            auto input = m_inputMeasurements(ctx);

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

            std::vector<Acts::SourceLink> sourceLinks;
            for (const auto& meas : input) {
                sourceLinks.push_back(meas.sourceLink);
            }

            auto gridLookup = m_cfg.sourceLinkGridConstructor->constructGrid(ctx.geoContext, sourceLinks);

            std::vector<typename Acts::Experimental::SeedTree> pathSeeds = 
                m_cfg.seeder->getSeeds(ctx.geoContext, gridLookup);

            Seeds outSeeds;
            auto me = 0.511 * Acts::UnitConstants::MeV;
            for (const auto& seed : pathSeeds) {
                Acts::Vector4 mPos4 = {seed.ipVertex.x(), seed.ipVertex.y(), seed.ipVertex.z(), 0};

                Acts::ActsScalar p = seed.ipP; 
                Acts::ActsScalar theta = std::acos(seed.ipDir.z()/p);
                Acts::ActsScalar phi = std::atan2(seed.ipDir.y(), seed.ipDir.x());

                Acts::CurvilinearTrackParameters ipParameters(
                    mPos4, phi, theta,
                    1_e / p, ipCov, 
                    Acts::ParticleHypothesis::electron());

                // Traverse the seed tree and collect the source links
                TryAllTracks tryAll;
                tryAll.tryAllTracks(seed.root, {});

                for (const auto& track : tryAll.tracks) {
                    outSeeds.push_back(Seed{
                        track,
                        ipParameters,
                        outSeeds.size()});
                }
            }

            std::cout << "PATH SEEDS SIZE: " << outSeeds.size() << std::endl;

            m_outputSeeds(ctx, std::move(outSeeds));

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
