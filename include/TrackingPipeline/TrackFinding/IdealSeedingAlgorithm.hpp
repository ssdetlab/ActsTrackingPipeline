#pragma once

#include "Acts/EventData/SourceLink.hpp"

#include <cstddef>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"

class IdealSeedingAlgorithm : public IAlgorithm {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Input source links
    std::string inputSourceLinks;
    /// Input sim clusters
    std::string inputSimClusters;
    /// Output seeds
    std::string outputSeeds;
    /// Lower cutoff on the seed size
    std::size_t minSeedSize;
    /// Higher cutoff on the seed size
    std::size_t maxSeedSize;
    /// Lower cutoff on the number of layers in a seed
    std::size_t minLayers;
    /// Higher cutoff on the number of layers in a seed
    std::size_t maxLayers;
  };

  /// @brief Constructor
  IdealSeedingAlgorithm(const Config& config, Acts::Logging::Level level);

  /// @brief Execute method
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// @brief Readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;

  ReadDataHandle<SimClusters> m_inputSimClusters{this, "InputSimClusters"};

  ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
      this, "InputSourceLinks"};

  WriteDataHandle<Seeds> m_outputSeeds{this, "OutputSeeds"};
};
