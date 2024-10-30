#pragma once

#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Simulation/Generators.hpp"
#include "TrackingPipeline/Clustering/IClusterFilter.hpp"

using namespace Acts::UnitLiterals;

class NoiseEmbeddingAlgorithm : public IAlgorithm {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// Noise generator
            std::shared_ptr<INoiseGenerator> noiseGenerator;
            /// Cluster Filter
            std::shared_ptr<IClusterFilter> clusterFilter;
            /// Detector geometry
            const Acts::Experimental::Detector* detector = nullptr;
            /// Input source links
            std::string inputSourceLinks = "Measurements";
            /// Input sim clusters
            std::string inputSimClusters = "SimClusters";
            /// Output source links
            std::string outputSourceLinks = "Measurements";
            /// Output sim clusters
            std::string outputSimClusters = "SimClusters";
        };

        /// @brief Constructor
        NoiseEmbeddingAlgorithm(Config config, Acts::Logging::Level level)
            : IAlgorithm("NoiseEmbeddingAlgorithm", level),
            m_cfg(std::move(config)) {
                if (m_cfg.detector == nullptr) {
                    throw std::invalid_argument("Detector is not set");
                }
                m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
                m_inputSimClusters.initialize(m_cfg.inputSimClusters);

                m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
                m_outputSimClusters.initialize(m_cfg.outputSimClusters);
        }
        ~NoiseEmbeddingAlgorithm() = default;

        /// @brief The execute method        
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            // Get the inputs from the context
            auto inputSourceLinks = m_inputSourceLinks(ctx);
            auto inputSimClusters = m_inputSimClusters(ctx);

            std::random_device rd;
            std::mt19937 gen(rd());

            // Dummy IP parameters
            Acts::CurvilinearTrackParameters ipParameters(
                Acts::Vector4::Zero(),
                Acts::Vector3::Zero(),
                1_e / 1_GeV,
                Acts::BoundSquareMatrix::Identity(),
                Acts::ParticleHypothesis::electron());

            int id = 0;
            std::vector<const Acts::Surface*> surfs;
            for (auto vol : m_cfg.detector->volumes()) {
                for (auto surf : vol->surfaces()) {
                    surfs.push_back(surf);
                }
            }

            auto noiseHits = m_cfg.noiseGenerator->gen(ctx.geoContext, gen, surfs);

            for (auto& hit : noiseHits) {
                hit.setEventId(ctx.eventNumber);
                hit.setIndex(inputSourceLinks.size());

                SimHit sm{
                    Acts::SourceLink(hit), 
                    Acts::BoundVector::Zero(), 
                    ipParameters,
                    -1,
                    id,
                    id};
                id--;
                SimCluster sc{
                    hit,
                    {sm},
                    false,
                    hit.index()};
                
                if (!m_cfg.clusterFilter->operator()(ctx.geoContext, sc)) {
                    continue;
                }

                inputSourceLinks.push_back(Acts::SourceLink(hit));
                inputSimClusters.push_back(sc);
            }

            m_outputSourceLinks(ctx, std::move(inputSourceLinks));
            m_outputSimClusters(ctx, std::move(inputSimClusters));

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;

        ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
            this, "InputSourceLinks"};

        ReadDataHandle<SimClusters> m_inputSimClusters{
            this, "InputSimClusters"};

        WriteDataHandle<std::vector<Acts::SourceLink>> m_outputSourceLinks{
            this, "OutputSourceLinks"};

        WriteDataHandle<SimClusters> m_outputSimClusters{
            this, "OutputSimClusters"};
};
