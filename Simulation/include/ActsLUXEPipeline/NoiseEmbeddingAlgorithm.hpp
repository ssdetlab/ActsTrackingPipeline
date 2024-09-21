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
                m_inputMeasurement.initialize(m_cfg.inputCollection);
                m_outputMeasurements.initialize(m_cfg.outputCollection);
        }
        ~NoiseEmbeddingAlgorithm() = default;

        /// @brief The execute method        
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            // Get the input seeds
            // from the context
            
            // auto start = std::chrono::system_clock::now();

            auto input = m_inputMeasurement(ctx);
            // SimMeasurements input;

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
                        SimMeasurement sm{
                            hit, 
                            Acts::BoundVector::Zero(),  
                            ipParameters,
                            id--};
                        input.push_back(sm);
                    }
                }
            }

            // auto end = std::chrono::system_clock::now();

            // std::cout << "Noise embedding took "
                // << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                // << "ms" << std::endl;

            m_outputMeasurements(ctx, std::move(input));

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;

        ReadDataHandle<SimMeasurements> m_inputMeasurement{
            this, "InputMeasurement"};

        WriteDataHandle<SimMeasurements> m_outputMeasurements{
            this, "OutputMeasurements"};
};
