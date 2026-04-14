#pragma once

#include <Acts/Definitions/Algebra.hpp>

#include <cstddef>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/ITrackParametersEstimator.hpp"

class TryAllSeedingAlgorithm : public IAlgorithm {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Input detector source links
    std::string inputSourceLinks;
    /// Input detector source links indices
    std::string inputSourceLinkIndices;
    /// Output seeds
    std::string outputSeeds;
    /// Output track parameters
    std::string outputTrackParameters;
    /// Track parameters estimator
    std::shared_ptr<ITrackParametersEstimator> trackParametersEstimator;
    /// Low cutoff on the number of layers in a seed
    std::size_t minLayers;
    /// High cutoff on the number of layers in a seed
    std::size_t maxLayers;
  };

  /// @brief Constructor
  TryAllSeedingAlgorithm(Config config, Acts::Logging::Level level);

  /// @brief The execute method
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Get readonly access to the config parameters
  const Config& config() const;

 private:
  /// Configuration
  Config m_cfg;

  /// Read data handles
  ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
      this, "InputSourceLinks"};

  ReadDataHandle<std::vector<std::size_t>> m_inputSourceLinkIndices{
      this, "InputSourceLinkIndices"};

  /// Write data handles
  WriteDataHandle<IndexSeeds> m_outputSeeds{this, "OutputSeeds"};

  WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_outputTrackParameters{this, "OutputTrackParameters"};
};
