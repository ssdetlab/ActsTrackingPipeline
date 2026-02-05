#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>

#include "TChain.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TVector3.h"
#include "TVectorD.h"
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
    std::string outputSeedsGuess;
    /// Output fitted seeds
    std::string outputSeedsEst;
    /// Output tracks (for cleaning)
    std::string outputTracks; 
    /// The names of the input files
    std::vector<std::string> filePaths;
    /// Name of the input tree
    std::string treeName;
    /// Chi2 cut
    double minChi2;
    double maxChi2;
    /// Merge into a single event flag
    bool mergeIntoOneEvent;
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

  /// WriteDataHandle for the guess seed data
  WriteDataHandle<Seeds> m_outputSeedsGuess{this, "SeedsGuess"};

  /// WriteDataHandle for the estimated seed data
  WriteDataHandle<Seeds> m_outputSeedsEst{this, "SeedsEst"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// Mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// Vector of {eventNr, entryMin, entryMax}
  std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> m_eventMap;

  /// The input tree name
  // TChain* m_chain = nullptr;
  TTree* m_chain = nullptr;
  TFile* m_file = nullptr;

 protected:
  /// Measurement hits
  std::vector<TVector3>* m_trackHitsGlobal = nullptr;
  std::vector<TVector2>* m_trackHitsLocal = nullptr;

  /// Covariances of the track hits
  std::vector<TMatrixD>* m_trackHitCovs = nullptr;

  /// Geometry ids of the track hits
  std::vector<std::size_t>* m_geometryIds = nullptr;

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
  std::size_t m_ndf;

  /// TrackId
  std::size_t m_trackId;

  /// EventId
  std::size_t m_eventId;

  /// PDG ID
  int m_pdgId;

  /// Charge
  int m_charge;

  /// Guessed bound track parameters
  TVectorD* m_boundTrackParametersGuess = nullptr;
  TMatrixD* m_boundTrackCovGuess = nullptr;

  /// KF predicted bound track parameters
  TVectorD* m_boundTrackParametersEst = nullptr;
  TMatrixD* m_boundTrackCovEst = nullptr;

  /// Initial guess of the momentum at the IP
  TLorentzVector* m_ipMomentumGuess = nullptr;

  /// Initial guess of the vertex at the IP
  TVector3* m_vertexGuess = nullptr;

  /// KF predicted momentum at the IP
  TLorentzVector* m_ipMomentumEst = nullptr;

  /// KF predicted vertex at the IP
  TVector3* m_vertexEst = nullptr;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
