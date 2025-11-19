#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <RtypesCore.h>

#include "DetectorEvent.hpp"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TrackingPipeline/Geometry/ApollonGeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

namespace ApollonIo {

class ApollonEventFilter;
class ApollonClusterFilter;

/// @brief ROOT file reader designed for the EUDAQ2 format
///
/// @note The events are assumed to be ordered
class ApollonRootDataReader : public IReader {
 public:
  using Hit = std::pair<std::size_t, std::size_t>;

  /// @brief The nested configuration struct
  struct Config {
    /// Collection with the measurement data
    std::string outputSourceLinks;
    /// The names of the input files
    std::vector<std::string> filePaths;
    /// Name of the input tree
    std::string treeName;
    /// The keys we have in the ROOT file
    std::string eventKey;
    /// Geometry ID scope
    int minGeoId;
    int maxGeoId;
    /// Surface map for high-precision local to global conversion
    std::map<Acts::GeometryIdentifier, const Acts::Surface*> surfaceMap;
  };

  ApollonRootDataReader(const ApollonRootDataReader&) = delete;
  ApollonRootDataReader(const ApollonRootDataReader&&) = delete;

  /// Constructor
  /// @param config The Configuration struct
  /// @param level The log level
  ApollonRootDataReader(const Config& config, Acts::Logging::Level level);

  /// Reader name() method
  virtual std::string name() const override { return "ApollonRootDataReader"; }

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

  /// WriteDataHandle for the observable data
  WriteDataHandle<std::vector<Acts::SourceLink>> m_outputSourceLinks{
      this, "OutputData"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// Mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// Vector of {eventNr, entryMin, entryMax}
  std::vector<std::tuple<uint32_t, std::size_t, std::size_t>> m_eventMap;

  /// The input tree name
  // TChain* m_chain = nullptr;
  TTree* m_chain = nullptr;
  TFile* m_file = nullptr;

  ApollonGeometry::GeometryOptions m_gOpt;

 protected:
  /// Detector event handle
  DetectorEvent* m_detEvent = nullptr;
  /// Event number handle
  ULong64_t m_eventId;
};

}  // namespace ApollonIo
