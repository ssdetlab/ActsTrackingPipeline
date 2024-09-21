#pragma once

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"

#include "ActsLUXEPipeline/E320SourceLinkGrid.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Seeding/PathSeeder.hpp"

class PathSeedingAlgorithm : public IAlgorithm {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// Path seeder
            std::shared_ptr<Acts::Experimental::PathSeeder> seeder;
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

            // auto start = std::chrono::system_clock::now();

            auto gridLookup = m_cfg.sourceLinkGridConstructor->constructGrid(ctx.geoContext, sourceLinks);

            std::vector<Acts::Experimental::PathSeeder::Seed> pathSeeds;
            m_cfg.seeder->getSeeds(ctx.geoContext, gridLookup, pathSeeds);

            // auto end = std::chrono::system_clock::now();

            // std::cout << "PATH SEEDING: Path seeding took "
                // << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                // << "ms" << std::endl;

            // start = std::chrono::system_clock::now();
            
            Seeds outSeeds;
            auto me = 0.511 * Acts::UnitConstants::MeV;
            for (int i = 0; i < pathSeeds.size(); i++) {
                const auto& seed = pathSeeds.at(i);
                Acts::Vector4 mPos4 = {seed.ipVertex.x(), seed.ipVertex.y(), seed.ipVertex.z(), 0};

                Acts::ActsScalar p = seed.ipP; 
                Acts::ActsScalar theta = std::acos(seed.ipDir.z()/p);
                Acts::ActsScalar phi = std::atan2(seed.ipDir.y(), seed.ipDir.x());

                Acts::CurvilinearTrackParameters ipParameters(
                    mPos4, phi, theta,
                    1_e / p, ipCov, 
                    Acts::ParticleHypothesis::electron());

                auto pivotSl = seed.sourceLinks.at(0);
                auto pivotSsl = pivotSl.get<SimpleSourceLink>();

                int trackId = 0;
                for (auto meas : input) {
                    auto sl = meas.sourceLink;
                    auto ssl = sl.get<SimpleSourceLink>();

                    if (pivotSsl == ssl) {
                        trackId = meas.trackId;
                        break;
                    }
                }

                outSeeds.push_back(Seed{
                    seed.sourceLinks,
                    ipParameters,
                    trackId});
            }

            // std::cout << "OUT SEEDS SIZE: " << outSeeds.size() << std::endl;

            m_outputSeeds(ctx, std::move(outSeeds));
            
            // end = std::chrono::system_clock::now();

            // std::cout << "PATH SEEDING: Conversion took "
                // << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                // << "ms" << std::endl;

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
