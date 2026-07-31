#pragma once

#include <cstddef>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Simulation/IMeasurementGenerator.hpp"

/// @brief Algorithm performing fast simulation of a particle
/// propagation and creating measurements on the geometry surfaces according to
/// the measurement generator implementation
class MeasurementsEmbeddingAlgorithm : public IAlgorithm {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Measurement generator
    std::shared_ptr<IMeasurementGenerator> measurementGenerator;
    /// Random number generator
    std::shared_ptr<RandomNumbers> randomNumberSvc;
    /// Input source links
    std::string inputSourceLinks;
    /// Input sim clusters
    std::string inputSimClusters;
    /// Input source links indices
    std::string inputSourceLinkIndices;
    /// Output source links
    std::string outputSourceLinks;
    /// Output sim clusters
    std::string outputSimClusters;
    /// Output source links indices
    std::string outputSourceLinkIndices;
    /// Number of measurements
    std::size_t nMeasurements;
  };

  /// @brief Constructor
  ///
  /// @param config configuration struct
  /// @param level logging level
  MeasurementsEmbeddingAlgorithm(const Config& config,
                                 Acts::Logging::Level level);

  /// @brief The execute method
  ///
  /// @param ctx current algorithm context
  ///
  /// @return algorithm process code
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// @brief Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;

  /// Read data handles
  ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
      this, "InputSourceLinks"};

  ReadDataHandle<SimClusters> m_inputSimClusters{this, "InputSimClusters"};

  ReadDataHandle<std::vector<std::size_t>> m_inputSourceLinksIndices{
      this, "InputSourceLinkIndices"};

  /// Write data handles
  WriteDataHandle<std::vector<Acts::SourceLink>> m_outputSourceLinks{
      this, "OutputSourceLinks"};

  WriteDataHandle<SimClusters> m_outputSimClusters{this, "OutputSimClusters"};

  WriteDataHandle<std::vector<std::size_t>> m_outputSourceLinkIndices{
      this, "OutputSourceLinkIndices"};
};
