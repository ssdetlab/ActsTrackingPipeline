#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <cstddef>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"

class E320SeedingAlgorithm : public IAlgorithm {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// HT seeder
    std::shared_ptr<HoughTransformSeeder> htSeeder;
    /// HT seeder options
    HoughTransformSeeder::Options htOptions;
    /// Input source links
    std::string inputSourceLinks;
    /// Output seeds
    std::string outputSeeds;
    /// Lower cutoff on the number of layers in a seed
    std::size_t minLayers;
    /// Higher cutoff on the number of layers in a seed
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

  Acts::Vector3 m_detFirstLayerPoint;
  Acts::Vector3 m_detFirstLayerNormal;
  Acts::Vector3 m_backShift;

  double m_dipoleLength;
  double m_dipoleFieldStrength;

  double m_xCorrectorLength;
  double m_xCorrectorFieldStrength;

  Acts::BoundSquareMatrix m_ipCov;

  ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
      this, "InputSourceLinks"};

  WriteDataHandle<Seeds> m_outputSeeds{this, "OutputSeeds"};
};
