#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Seeding/PathSeeder.hpp"

#include <cstddef>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/E320SourceLinkGridConstructor.hpp"

/// @brief Algorithm performing path seeding
///
/// Algorithm encapsulates PathSeeder to provide
/// seed estimates to the subsequent track estimation algorithms
class PathSeedingAlgorithm : public IAlgorithm {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Path seeder
    std::shared_ptr<Acts::PathSeeder> seeder;
    /// SourceLink grid
    std::shared_ptr<E320SourceLinkGridConstructor> sourceLinkGridConstructor;
    /// Input source links
    std::string inputSourceLinks = "SourceLink";
    /// Output seeds
    std::string outputSeeds = "Seed";
    /// Lower cutoff on the seed size
    std::size_t minSeedSize;
    /// Higher cutoff on the seed size
    std::size_t maxSeedSize;
  };

  /// @brief Constructor
  PathSeedingAlgorithm(const Config& config, Acts::Logging::Level level);

  /// @brief Execute method
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// @brief Readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;

  ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
      this, "InputSourceLinks"};

  WriteDataHandle<Seeds> m_outputSeeds{this, "OutputSeeds"};
};
