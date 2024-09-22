#pragma once

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/ProxyAccessor.hpp"

template <typename propagator_t,
typename trajectory_t = Acts::VectorMultiTrajectory,
typename container_t = Acts::VectorTrackContainer>
class TrackFittingAlgorithm : public IAlgorithm {
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
        TrackFittingAlgorithm(Config config, Acts::Logging::Level level)
            : IAlgorithm("TrackFittingAlgorithm", level),
            m_cfg(std::move(config)) {
                m_inputSeeds.initialize(m_cfg.inputCollection);
                m_outputTracks.initialize(m_cfg.outputCollection);
        }
        ~TrackFittingAlgorithm() = default;

        /// @brief The execute method        
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            // auto start = std::chrono::system_clock::now();

            // Get the input seeds
            // from the context
            auto input = m_inputSeeds(ctx);

            auto trackContainer = std::make_shared<container_t>();
            auto trackStateContainer = std::make_shared<trajectory_t>();
            Acts::TrackContainer tracks(trackContainer, trackStateContainer);

            std::vector<std::int32_t> trackIds;

            for (const auto& seed : input) {
                auto start = seed.ipParameters;
                auto sourceLinks = seed.sourceLinks;

                auto res = m_cfg.fitter.fit(sourceLinks.begin(), sourceLinks.end(), 
                    start, m_cfg.kfOptions, tracks);
                
                if (!res.ok()) {
                    ACTS_ERROR("Track fitting failed");
                    continue;
                }

                trackIds.push_back(seed.trackId);
            }
            auto outTracks = Tracks<container_t, trajectory_t>{
                tracks, trackIds};

            m_outputTracks(ctx, std::move(outTracks));

            // auto end = std::chrono::system_clock::now();

            // std::cout << "Track fitting took "
                // << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                // << "ms" << std::endl;

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;

        ReadDataHandle<Seeds> m_inputSeeds
            {this, "InputSeeds"};

        WriteDataHandle<Tracks<
            container_t, trajectory_t>> m_outputTracks
            {this, "OutputTracks"};
};
