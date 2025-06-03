#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <set>

#include <RtypesCore.h>

#include "TChain.h"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Io/detail/AcceleratorState.hpp"
#include "TrackingPipeline/Io/detail/DetectorEvent.hpp"

namespace E320Io {

class E320EventFilter;
class E320ClusterFilter;

/// @brief ROOT file reader designed for the EUDAQ2 format
///
/// @note The events are assumed to be ordered
class E320RootDataReader : public IReader {
 public:
  using Hit = std::pair<std::size_t, std::size_t>;

  /// @brief The nested configuration struct
  struct Config {
    /// Cluster filter
    std::shared_ptr<E320ClusterFilter> clusterFilter = nullptr;
    /// Event filter
    std::shared_ptr<E320EventFilter> dataFilter = nullptr;
    

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

  E320RootDataReader(const E320RootDataReader&) = delete;
  E320RootDataReader(const E320RootDataReader&&) = delete;

  /// Constructor
  /// @param config The Configuration struct
  /// @param level The log level
  E320RootDataReader(const Config& config, Acts::Logging::Level level);

  /// Reader name() method
  virtual std::string name() const override { return "E320RootDataReader"; }

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
  DetectorEvent* m_detEvent = nullptr;
  /// Stable accelerator state
  AcceleratorState m_stableAcceleratorState;
  /// Event number handle
  ULong64_t m_eventId;
};

class E320EventFilter {
 public:
  struct Config {
    /// Maximum allowed turnaround time in ms
    double maxTurnaroundTime;
    /// Maximum allowed number of standard deviations
    /// from the mean of the accelerator state
    double nStdDevs;
  };

  E320EventFilter(const E320EventFilter&) = delete;
  E320EventFilter(const E320EventFilter&&) = delete;

  /// Constructor
  /// @param config The Configuration struct
  /// @param level The log level
  E320EventFilter(const Config& config, Acts::Logging::Level level);

  void appendAcceleratorState(const EpicsFrame& epicsFrame);

  void finalizeAcceleratorState();

  bool checkCuts(const DetectorEvent& detEvent);

 private:
  Config m_cfg;
  AcceleratorState m_state;
  std::size_t m_stateSize = 0;
};

class E320ClusterFilter {
 public:
  double a1 = 0.11672788;
  double b1 = 5.3475191;

  double a2 = 0.13025286;
  double b2 = 3.9499999;

  bool operator()(double hitX, double hitY) {
    if (hitY > a2 * hitX + b2 && hitY < a1 * hitX + b1) {
      return true;
    } else {
      return false;
    }
  }
};
}  // namespace E320Io
