#pragma once

#include <cstddef>

#include "TrackingPipeline/EventData/SimCluster.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

/// @brief dummy reader enabling Sequencer's event loop but
/// providing no data
class DummyReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Output source links
    std::string outputSourceLinks;
    /// Output sim clusters
    std::string outputSimClusters;
    /// Output source links indices
    std::string outputSourceLinkIndices;
    /// Number of events
    std::size_t nEvents;
  };

  /// @brief Constructor
  ///
  /// @param config configuration struct
  explicit DummyReader(const Config& config);

  /// @brief The execute method
  ///
  /// @param ctx current algorithm context
  ///
  /// @return algorithm process code
  ProcessCode read(const AlgorithmContext& ctx) override;

  /// @brief Provide range of available events
  std::pair<std::size_t, std::size_t> availableEvents() const override;

  /// @brief get reader's name
  std::string name() const override { return "DummyReader"; }

  /// @brief Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;

  /// Write data handles
  WriteDataHandle<std::vector<Acts::SourceLink>> m_outputSourceLinks{
      this, "OutputSourceLinks"};

  WriteDataHandle<SimClusters> m_outputSimClusters{this, "OutputSimClusters"};

  WriteDataHandle<std::vector<std::size_t>> m_outputSourceLinkIndices{
      this, "OutputSourceLinkIndices"};
};
