#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <RtypesCore.h>

#include "TChain.h"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Io/detail/DetectorEvent.hpp"

namespace E320Io {

class E320EventFilter;
class E320ClusterFilter;

/// @brief ROOT file reader designed for the EUDAQ2 format
///
/// @note The events are assumed to be ordered
class RootNoamSplitDataReader : public IReader {
 public:
  using Hit = std::pair<std::size_t, std::size_t>;

  /// @brief The nested configuration struct
  struct Config {
    /// Surface accessor
    Acts::SourceLinkSurfaceAccessor surfaceAccessor;


    /// Collection with the measurement data
    std::string outputSourceLinks;
    /// The names of the input files
    std::vector<std::string> filePaths;
    /// Name of the input tree
    std::string treeName = "MyTree";
    /// The keys we have in the ROOT file
    std::string eventKey = "event";
    /// Number of triggers to skip
    std::size_t skip = 0;
  };

  RootNoamSplitDataReader(const RootNoamSplitDataReader&) = delete;
  RootNoamSplitDataReader(const RootNoamSplitDataReader&&) = delete;

  /// Constructor
  /// @param config The Configuration struct
  /// @param level The log level
  RootNoamSplitDataReader(const Config& config, Acts::Logging::Level level);

  /// Reader name() method
  virtual std::string name() const override { return "RootNoamSplitDataReader"; }

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
  TChain* m_chain = nullptr;

  E320Geometry::GeometryOptions m_gOpt;

 protected:
  /// Detector event handle
  event* m_detEvent = nullptr;
};

}  // namespace E320Io

