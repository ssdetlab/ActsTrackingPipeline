#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"

#include <cstddef>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"

class E320SeedingAlgorithm : public IAlgorithm {
 public:
  enum struct PropagationDirection : int { forward = 0, backward = 1 };

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
    /// Reference surface
    const Acts::Surface* referenceSurface;
    /// Initial track state covariance prior
    Acts::BoundMatrix originCov;
    /// Lower cutoff on the number of layers in a seed
    std::size_t minLayers;
    /// Higher cutoff on the number of layers in a seed
    std::size_t maxLayers;
    /// BPM ids
    std::vector<int> bpmIds;
    /// Forward or backward propagation parameters
    PropagationDirection propDirection;
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

  Acts::BoundMatrix transportCovToReference(
      const Acts::GeometryContext& gctx, const Acts::Vector3& refSurfacePoint,
      const Acts::Vector3& point, const Acts::Vector3& dir,
      const Acts::ActsSquareMatrix<6>& cov) const;

  double m_dipoleLength;
  double m_dipoleFieldStrength;

  ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
      this, "InputSourceLinks"};

  WriteDataHandle<Seeds> m_outputSeeds{this, "OutputSeeds"};
};
