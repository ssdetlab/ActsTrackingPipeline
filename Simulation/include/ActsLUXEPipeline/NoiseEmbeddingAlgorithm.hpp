#pragma once

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"
#include "ActsLUXEPipeline/Generators.hpp"

using namespace Acts::UnitLiterals;

class NoiseEmbeddingAlgorithm : public IAlgorithm {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// Noise generator
            std::shared_ptr<INoiseGenerator> noiseGenerator;
            /// Detector geometry
            const Acts::Experimental::Detector* detector = nullptr;
            /// The input collection
            std::string inputCollection = "Measurements";
            /// The output collection
            std::string outputCollection = "Measurements";
        };

        /// @brief Constructor
        NoiseEmbeddingAlgorithm(Config config, Acts::Logging::Level level)
            : IAlgorithm("NoiseEmbeddingAlgorithm", level),
            m_cfg(std::move(config)) {
                if (m_cfg.detector == nullptr) {
                    throw std::invalid_argument("Detector is not set");
                }
                m_inputSourceLinks.initialize(m_cfg.inputCollection);
                m_outputSourceLinks.initialize(m_cfg.outputCollection);
        }
        ~NoiseEmbeddingAlgorithm() = default;

        /// @brief The execute method        
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            // Get the inputs from the context
            auto inputSourceLinks = m_inputSourceLinks(ctx);
            auto inputSimHits = m_inputSimClusters(ctx);

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
            for (auto vol : m_cfg.detector->volumes()) {
                for (auto surf : vol->surfaces()) {
                    auto noiseHits = m_cfg.noiseGenerator->gen(gen, surf);
                    for (auto& hit : noiseHits) {
                        hit.setIndex(inputSourceLinks.size());
                        inputSourceLinks.push_back(Acts::SourceLink(hit));

                        SimHit sm{
                            Acts::SourceLink(hit), 
                            Acts::BoundVector::Zero(), 
                            ipParameters,
                            -1,
                            id--,
                            id--};
                        inputSimHits.push_back(sm);
                    }
                }
            }

            m_outputSourceLinks(ctx, std::move(inputSourceLinks));
            m_outputSimClusters(ctx, std::move(inputSimHits));

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;

        ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
            this, "InputSourceLinks"};

        ReadDataHandle<SimHits> m_inputSimClusters{
            this, "InputSimClusters"};

        WriteDataHandle<std::vector<Acts::SourceLink>> m_outputSourceLinks{
            this, "OutputSourceLinks"};

        WriteDataHandle<SimHits> m_outputSimClusters{
            this, "OutputSimClusters"};
};
