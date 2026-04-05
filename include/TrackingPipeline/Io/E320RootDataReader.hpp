#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <unordered_map>

#include <RtypesCore.h>

#include "DetectorEvent.hpp"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

namespace E320 {

/// @brief ROOT file reader designed for the EUDAQ2 format
class E320RootDataReader : public IReader {
 public:
  using Hit = std::pair<std::size_t, std::size_t>;

  /// @brief The nested configuration struct
  struct Config {
    /// Collection with the detector measurement data
    std::string outputSourceLinks;
    /// Collection with the detector measurement data indices
    std::string outputDetSourceLinkIndices;
    /// Collection with the BPM measurement data indices
    std::string outputBpmSourceLinkIndices;
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
    std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
        surfaceMap;
  };

  E320RootDataReader(const E320RootDataReader&) = delete;
  E320RootDataReader(const E320RootDataReader&&) = delete;

  /// @brief constructor
  E320RootDataReader(const Config& config, Acts::Logging::Level level);

  /// @brief get reader name
  std::string name() const override { return "E320RootDataReader"; }

  /// @brief return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const final;

  /// @brief read out data from the input stream
  ProcessCode read(const AlgorithmContext& ctx) override;

  /// @brief readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// @brief private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  /// WriteDataHandle for the observable data
  WriteDataHandle<std::vector<Acts::SourceLink>> m_outputSourceLinks{
      this, "OutputSourceLinks"};

  WriteDataHandle<std::vector<std::size_t>> m_outputDetSourceLinksIndices{
      this, "OutputDetSourceLinksIndices"};

  WriteDataHandle<std::vector<std::size_t>> m_outputBpmSourceLinksIndices{
      this, "OutputBpmSourceLinksIndices"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// Mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// Vector of {eventNr, entryMin, entryMax}
  std::vector<std::tuple<uint32_t, std::size_t, std::size_t>> m_eventMap;

  std::unordered_map<std::uint8_t, std::size_t> m_geoIdMap{
      {0, 18}, {2, 16},    {4, 14},    {6, 12},
      {8, 10}, {3218, 41}, {3265, 42}, {3315, 43}};

  /// The input tree name
  TTree* m_tree = nullptr;
  TFile* m_file = nullptr;
  TChain* m_chainOwner = nullptr;

  E320::GeometryOptions m_gOpt;

 protected:
  /// Detector event handle
  E320Io::DetectorEvent* m_detEvent = nullptr;
  /// Event number handle
  ULong64_t m_eventId = 0;
};

}  // namespace E320
