#pragma once

#include "Acts/EventData/SourceLink.hpp"

#include <cstddef>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"
#include "TrackingPipeline/TrackFinding/ITrackParametersEstimator.hpp"

namespace E320 {

class E320SeedingAlgorithm : public IAlgorithm {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// HT seeder
    std::shared_ptr<HoughTransformSeeder> htSeeder;
    /// HT seeder options
    HoughTransformSeeder::Options htOptions;
    /// Input detector source links
    std::string inputSourceLinks;
    /// Input detector source links indices
    std::string inputDetSourceLinkIndices;
    /// Input BPM source links indices
    std::string inputBpmSourceLinkIndices;
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
  E320SeedingAlgorithm(const Config& config, Acts::Logging::Level level);

  /// @brief Execute method
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// @brief Readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;

  /// Read data handles
  ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
      this, "InputSourceLinks"};

  ReadDataHandle<std::vector<std::size_t>> m_inputDetSourceLinkIndices{
      this, "InputDetSourceLinkIndices"};

  ReadDataHandle<std::vector<std::size_t>> m_inputBpmSourceLinkIndices{
      this, "InputBpmSourceLinkIndices"};

  /// Write data handles
  WriteDataHandle<IndexSeeds> m_outputSeeds{this, "OutputSeeds"};

  WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_outputTrackParameters{this, "OutputTrackParameters"};
};

}  // namespace E320
