#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/ProxyAccessor.hpp"

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"
#include "ActsLUXEPipeline/LUXEDataContainers.hpp"

template <typename propagator_t, 
typename trajectory_t = Acts::VectorMultiTrajectory>
class TrackFitter : public IAlgorithm {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// The input collection
            std::string inputCollection = "SourceLink";
            /// The output collection
            std::string outputCollection = "Track";
            /// KF fitter
            const Acts::KalmanFitter<propagator_t, trajectory_t>& fitter;
            /// KF options
            Acts::KalmanFitterOptions<trajectory_t> kfOptions;
        };

        /// @brief Constructor
        TrackFitter(Config config, Acts::Logging::Level level)
            : IAlgorithm("TrackFitter", level),
            m_cfg(std::move(config)) {
                m_inputSeeds.initialize(m_cfg.inputCollection);
                m_outputTracks.initialize(m_cfg.outputCollection);
        }
        ~TrackFitter() = default;

        /// @brief The execute method        
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            // Get the input seeds
            // from the context
            auto input = m_inputSeeds(ctx);

            Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                                        trajectory_t{}};

            for (const auto& seed : input) {
                auto start = seed.ipParameters;
                auto sourceLinks = seed.sourceLinks;
                auto res = m_cfg.fitter.fit(sourceLinks.begin(), sourceLinks.end(), 
                    start, m_cfg.kfOptions, tracks);

                break;
            }

            // m_outputTracks(ctx, std::move(res.value()));

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;

        ReadDataHandle<LUXEDataContainer::Seeds> m_inputSeeds
            {this, "InputSeeds"};

        WriteDataHandle<LUXEDataContainer::Seeds> m_outputTracks
            {this, "OutputTracks"};
};
