#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>

#include "TChain.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TrackingPipeline/Clustering/IClusterFilter.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

/// @brief Intermediate generalization of the
/// ROOT file reader to be inhereted from by the
/// readers for the specific tree structures,
/// data types and geometries
///
/// @note The events are assumed to be ordered
class RootSimDataReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Cluster filter
    std::shared_ptr<IClusterFilter> clusterFilter = nullptr;
    /// Collection with the measurement data
    std::string outputSourceLinks;
    /// Collection with the sim clusters data
    std::string outputSimClusters;
    /// Name of the input tree
    std::string treeName = "clusters";
    /// The names of the input files
    std::vector<std::string> filePaths;
    /// The keys we have in the ROOT file
    std::vector<const char*> vVector3Keys;
    std::vector<const char*> vVector2Keys;
    std::vector<const char*> vector3Keys;
    std::vector<const char*> vLorentzKeys;
    std::vector<const char*> vIntKeys;
    std::vector<const char*> vDoubleKeys;
    std::vector<const char*> intKeys;
  };

  RootSimDataReader(const RootSimDataReader&) = delete;
  RootSimDataReader(const RootSimDataReader&&) = delete;

  /// Constructor
  /// @param config The Configuration struct
  /// @param level The log level
  RootSimDataReader(const Config& config, Acts::Logging::Level level);

  /// Reader name() method
  virtual std::string name() const override { return "RootSimDataReader"; }

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const final;

  /// Read out data from the input stream
  ProcessCode read(const AlgorithmContext& context) override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  /// Prepare the measurements
  virtual inline void prepareMeasurements(
      const AlgorithmContext& context,
      std::vector<Acts::SourceLink>* sourceLinks,
      SimClusters* clusters) const = 0;

  /// WriteDataHandle for the observable data
  WriteDataHandle<std::vector<Acts::SourceLink>> m_outputSourceLinks{
      this, "OutputObsData"};

  /// WriteDataHandle for the simulation data
  WriteDataHandle<SimClusters> m_outputSimClusters{this, "OutputSimClusters"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// Vector of {eventNr, entryMin, entryMax}
  std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> m_eventMap;

  /// The input tree name
  TChain* m_chain = nullptr;

 protected:
  /// The exausitive list of columns
  std::unordered_map<std::string_view, int> m_intColumns;
  std::unordered_map<std::string_view, std::vector<int>*> m_vIntColumns;

  std::unordered_map<std::string_view, std::double_t> m_doubleColumns;
  std::unordered_map<std::string_view, std::vector<std::double_t>*>
      m_vDoubleColumns;

  std::unordered_map<std::string_view, TVector2*> m_vector2Columns;
  std::unordered_map<std::string_view, std::vector<TVector2>*>
      m_vVector2Columns;

  std::unordered_map<std::string_view, TVector3*> m_vector3Columns;
  std::unordered_map<std::string_view, std::vector<TVector3>*>
      m_vVector3Columns;

  std::unordered_map<std::string_view, TLorentzVector*> m_lorentzColumns;
  std::unordered_map<std::string_view, std::vector<TLorentzVector>*>
      m_vLorentzColumns;
};
