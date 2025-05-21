#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>

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

class RootTrackReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Output source links
    std::string outputMeasurements;
    /// Output seeds
    std::string outputSeeds;
    /// The names of the input files
    std::vector<std::string> filePaths;
    /// Name of the input tree
    std::string treeName;
    /// Batch flag
    bool batch;
    /// Batch size for the tracks
    std::size_t batchSize;
    /// Covariance annealing factor
    double covAnnealingFactor;
  };

  RootTrackReader(const RootTrackReader&) = delete;
  RootTrackReader(const RootTrackReader&&) = delete;

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  RootTrackReader(const Config& config, Acts::Logging::Level level);

  /// Writer name() method
  std::string name() const override { return "RootTrackReader"; }

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

  /// WriteDataHandle for the observable data
  WriteDataHandle<std::vector<Acts::SourceLink>> m_outputSourceLinks{
      this, "OutputData"};

  /// WriteDataHandle for the sim cluster data
  WriteDataHandle<Seeds> m_outputSeeds{this, "Seeds"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// Mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// Vector of {eventNr, entryMin, entryMax}
  std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> m_eventMap;

  /// The input tree name
  TChain* m_chain = nullptr;

 protected:
  /// Measurement hits
  std::vector<TVector3>* m_trackHitsGlobal = nullptr;
  std::vector<TVector2>* m_trackHitsLocal = nullptr;

  /// Covariances of the track hits
  std::vector<TMatrixD>* m_trackHitCovs = nullptr;

  /// Geometry ids of the track hits
  std::vector<int>* m_geometryIds = nullptr;

  /// KF predicted track hits
  std::vector<TVector3>* m_predictedTrackHitsGlobal = nullptr;
  std::vector<TVector3>* m_filteredTrackHitsGlobal = nullptr;
  std::vector<TVector3>* m_smoothedTrackHitsGlobal = nullptr;

  std::vector<TVector2>* m_predictedTrackHitsLocal = nullptr;
  std::vector<TVector2>* m_filteredTrackHitsLocal = nullptr;
  std::vector<TVector2>* m_smoothedTrackHitsLocal = nullptr;

  /// KF residuals with respect to the measurements
  std::vector<TVector2>* m_predictedResiduals = nullptr;
  std::vector<TVector2>* m_filteredResiduals = nullptr;
  std::vector<TVector2>* m_smoothedResiduals = nullptr;

  /// KF pulls with respect to the measurements
  std::vector<TVector2>* m_predictedPulls = nullptr;
  std::vector<TVector2>* m_filteredPulls = nullptr;
  std::vector<TVector2>* m_smoothedPulls = nullptr;

  /// Chi2 of the track
  /// with respect ot the
  /// measurement
  double m_chi2Predicted;
  double m_chi2Filtered;
  double m_chi2Smoothed;

  /// Number of degrees of freedom
  /// of the track
  int m_ndf;

  /// TrackId
  std::vector<int>* m_trackId = nullptr;

  /// EventId
  int m_eventId;

  /// Initial guess of the momentum at the IP
  TLorentzVector* m_ipMomentumGuess = nullptr;
  TVector3* m_vertexGuess = nullptr;

  /// KF predicted momentum at the IP
  TLorentzVector* m_ipMomentumEst = nullptr;
  TVector3* m_ipMomentumError = nullptr;
  TVector3* m_vertexEst = nullptr;
  TVector3* m_vertexError = nullptr;

  /// PDG ID
  int m_pdgId;

  /// Charge
  int m_charge;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};

