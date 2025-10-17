#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <cstddef>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"

class ApollonSeedingAlgorithm : public IAlgorithm {
 public:
  enum class SeedingScope : int {
    detector1,
    detector2,
    fullDetector
  };

  /// @brief The nested configuration struct
  struct Config {
    /// HT seeder
    std::shared_ptr<HoughTransformSeeder> htSeeder;
    /// Input source links
    std::string inputSourceLinks;
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
    /// Beginning of the energy scan range
    double minScanEnergy;
    /// End of the energy scan range
    double maxScanEnergy;
    /// Energy scan step
    double energyScanStep;
    /// Max distance at which the seed connection is accepted
    double maxConnectionDistance;
    /// Scope of the seeding
    SeedingScope scope;
  };

  /// @brief Constructor
  ApollonSeedingAlgorithm(const Config& config, Acts::Logging::Level level);

  /// @brief Execute method
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// @brief Readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;

  Acts::Vector3 m_det1FirstLayerPoint;
  Acts::Vector3 m_det1FirstLayerNormal;

  Acts::Vector3 m_det2FirstLayerPoint;
  Acts::Vector3 m_det2FirstLayerNormal;

  Acts::Vector3 m_det2LastLayerPoint;
  Acts::Vector3 m_det2LastLayerNormal;

  Acts::Vector3 m_dipoleEntrancePoint;
  Acts::Vector3 m_dipoleEntranceNormal;

  Acts::Vector3 m_dipoleExitPoint;
  Acts::Vector3 m_dipoleExitNormal;

  double m_dAngleMax;
  double m_dOrthoMax;

  Seeds scanEnergy(
      const std::vector<HoughTransformSeeder::HTSeed>& det1Seeds,
      const std::vector<HoughTransformSeeder::HTSeed>& det2Seeds) const;

  Acts::BoundSquareMatrix m_ipCov;

  ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
      this, "InputSourceLinks"};

  WriteDataHandle<Seeds> m_outputSeeds{this, "OutputSeeds"};
};
