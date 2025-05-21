#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "RtypesCore.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVector3.h"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

class RootSimClusterReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Output source links
    std::string outputMeasurements;
    /// Output sim clusters
    std::string outputSimClusters;
    /// The names of the input files
    std::vector<std::string> filePaths;
    /// Name of the input tree
    std::string treeName = "MyTree";
  };

  RootSimClusterReader(const RootSimClusterReader&) = delete;
  RootSimClusterReader(const RootSimClusterReader&&) = delete;

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  RootSimClusterReader(const Config& config, Acts::Logging::Level level);

  /// Writer name() method
  std::string name() const override { return "RootSimClusterReader"; }

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const final;

  /// Write out data to the input stream
  ProcessCode read(const AlgorithmContext& ctx) override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  /// WriteDataHandle for the sim data
  WriteDataHandle<SimClusters> m_outputSimClusters{this, "SimClusters"};

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

 protected:
  TVector2* m_geoCenterLocal = nullptr;
  TMatrixD* m_cov = nullptr;

  std::size_t m_geoId;
  std::size_t m_eventId;

  int m_isSignal;

  /// Measurement hits
  std::vector<TVector2>* m_trackHits = nullptr;

  std::vector<int>* m_trackId = nullptr;
  std::vector<int>* m_parentTrackId = nullptr;
  std::vector<int>* m_runId = nullptr;

  std::vector<TLorentzVector>* m_originMomentum = nullptr;
  std::vector<TLorentzVector>* m_onSurfaceMomentum = nullptr;
  std::vector<TVector3>* m_vertex = nullptr;

  std::vector<int>* m_charge = nullptr;
  std::vector<int>* m_pdgId = nullptr;
};
