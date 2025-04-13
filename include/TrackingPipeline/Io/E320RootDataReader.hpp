#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <RtypesCore.h>

#include "TChain.h"
#include "TrackingPipeline/Clustering/IClusterFilter.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Io/detail/EudaqDetectorEvent.h"

namespace E320Io {

/// @brief ROOT file reader designed for the EUDAQ2 format
///
/// @note The events are assumed to be ordered
class E320RootDataReader : public IReader {
 public:
  using Hit = std::pair<std::size_t, std::size_t>;

  /// @brief The nested configuration struct
  struct Config {
    /// Cluster filter
    std::shared_ptr<IClusterFilter> clusterFilter = nullptr;
    /// Collection with the measurement data
    std::string outputSourceLinks;
    /// Name of the input tree
    std::string treeName = "MyTree";
    /// The names of the input files
    std::vector<std::string> filePaths;
    /// The keys we have in the ROOT file
    std::string eventKey = "event";
  };

  struct Cluster {
    double xCenter;
    double yCenter;

    std::size_t dx;
    std::size_t dy;
    std::size_t cl_size;
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
      this, "OutputObsData"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// Mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// Vector of {eventNr, entryMin, entryMax}
  std::vector<std::tuple<uint32_t, std::size_t, std::size_t>> m_eventMap;

  /// The input tree name
  TChain* m_chain = nullptr;

  E320Geometry::GeometryOptions m_gOpt;

  std::vector<Cluster> getClusters(std::set<Hit>& hits);

  void bfsClustering(std::vector<Hit>& clusterPixels, Hit& pivot,
                     std::set<Hit>& pixels);

 protected:
  /// Detector event handle
  detector_event* m_detEvent = nullptr;
  /// Event number handle
  ULong64_t m_eventId;
};

}  // namespace E320Io
