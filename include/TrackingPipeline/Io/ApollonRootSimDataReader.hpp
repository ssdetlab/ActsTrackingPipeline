#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>

#include "TChain.h"
#include "TVector3.h"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

namespace ApollonIo {

class ApollonRootSimDataReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Collection with the measurement data
    std::string outputSourceLinks;
    /// Collection with the sim clusters data
    std::string outputSimClusters;
    /// Name of the input tree
    std::string treeName;
    /// The names of the input files
    std::vector<std::string> filePaths;
    /// Whether to split data on G4 event on run basis
    bool eventSplit;
  };

  ApollonRootSimDataReader(const ApollonRootSimDataReader&) = delete;
  ApollonRootSimDataReader(const ApollonRootSimDataReader&&) = delete;

  /// Constructor
  /// @param config The Configuration struct
  /// @param level The log level
  ApollonRootSimDataReader(const Config& config, Acts::Logging::Level level);

  /// Reader name() method
  virtual std::string name() const override {
    return "ApollonRootSimDataReader";
  }

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

  Acts::BoundSquareMatrix m_ipCov;
  Acts::SquareMatrix2 m_hitCov;
  Acts::RotationMatrix3 m_setupRotation;
  Acts::Vector3 m_ipCorrection;

  /// WriteDataHandle for the observable data
  WriteDataHandle<std::vector<Acts::SourceLink>> m_outputSourceLinks{
      this, "OutputSourceLinks"};

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
  std::vector<int>* m_geoIdVal = nullptr;
  std::vector<int>* m_isSignalFlag = nullptr;

  std::vector<int>* m_trackId = nullptr;
  std::vector<int>* m_parentTrackId = nullptr;
  int m_eventId;
  int m_runId;

  std::vector<TVector3>* m_hitPosGlobal = nullptr;
  std::vector<TVector2>* m_hitPosLocal = nullptr;

  std::vector<TVector3>* m_hitMomDir = nullptr;
  std::vector<double>* m_hitE = nullptr;
  std::vector<double>* m_hitP = nullptr;
  std::vector<double>* m_eDep = nullptr;

  std::vector<TVector3>* m_ipMomDir = nullptr;
  std::vector<double>* m_ipE = nullptr;
  std::vector<double>* m_ipP = nullptr;
  std::vector<TVector3>* m_vertices = nullptr;
};

}  // namespace ApollonIo
