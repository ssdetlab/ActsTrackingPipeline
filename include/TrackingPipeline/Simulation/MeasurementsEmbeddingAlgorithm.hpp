#pragma once

#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Simulation/IMeasurementGenerator.hpp"
#include "TrackingPipeline/Clustering/IClusterFilter.hpp"

#include <cstddef>

using namespace Acts::UnitLiterals;

class MeasurementsEmbeddingAlgorithm : public IAlgorithm {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// Measurement generator
            std::shared_ptr<IMeasurementGenerator> measurementGenerator;
            /// Cluster Filter
            std::shared_ptr<IClusterFilter> clusterFilter = nullptr;
            /// Random number generator
            std::shared_ptr<RandomNumbers> randomNumberSvc;
            /// Input source links
            std::string inputSourceLinks;
            /// Input sim clusters
            std::string inputSimClusters;
            /// Output source links
            std::string outputSourceLinks;
            /// Output sim clusters
            std::string outputSimClusters;
            /// Number of measurements
            std::size_t nMeasurements;
        };

        /// @brief Constructor
        MeasurementsEmbeddingAlgorithm(
            const Config& config, Acts::Logging::Level level);

        /// @brief The execute method        
        ProcessCode execute(const AlgorithmContext& ctx) const override; 

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
