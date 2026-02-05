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

class RootSimTrackReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Output source links
    std::string outputMeasurements;
    /// Output sim clusters
    std::string outputSimClusters;
    /// Output seeds
    std::string outputSeedsGuess;
    /// Output fitted seeds
    std::string outputSeedsEst;
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

  RootSimTrackReader(const RootSimTrackReader&) = delete;
  RootSimTrackReader(const RootSimTrackReader&&) = delete;

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  RootSimTrackReader(const Config& config, Acts::Logging::Level level);

  /// Writer name() method
  std::string name() const override { return "RootSimTrackReader"; }

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
  WriteDataHandle<SimClusters> m_outputSimClusters{this, "SimClusters"};

  /// WriteDataHandle for the seed data
  WriteDataHandle<Seeds> m_outputSeedsGuess{this, "SeedsGuess"};

  /// WriteDataHandle for the fitted seed data
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
  /// True hits
  std::vector<TVector3>* m_trueTrackHitsGlobal = nullptr;
  std::vector<TVector2>* m_trueTrackHitsLocal = nullptr;
  std::vector<TLorentzVector>* m_onSurfaceMomentumTruth = nullptr;
  std::vector<int>* m_isSignal = nullptr;

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

  /// KF residuals with respect to the true hits
  std::vector<TVector2>* m_truePredictedResiduals = nullptr;
  std::vector<TVector2>* m_trueFilteredResiduals = nullptr;
  std::vector<TVector2>* m_trueSmoothedResiduals = nullptr;

  /// KF residuals with respect to the measurements
  std::vector<TVector2>* m_predictedResiduals = nullptr;
  std::vector<TVector2>* m_filteredResiduals = nullptr;
  std::vector<TVector2>* m_smoothedResiduals = nullptr;

  /// KF pulls with respect to the true hits
  std::vector<TVector2>* m_truePredictedPulls = nullptr;
  std::vector<TVector2>* m_trueFilteredPulls = nullptr;
  std::vector<TVector2>* m_trueSmoothedPulls = nullptr;

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
  std::vector<std::size_t>* m_stateTrackId = nullptr;
  std::vector<std::size_t>* m_stateParentTrackId = nullptr;
  std::vector<std::size_t>* m_stateRunId = nullptr;

  std::size_t m_trackId;
  std::size_t m_parentTrackId;
  std::size_t m_runId;

  /// EventId
  std::size_t m_eventId;

  /// True track size
  std::size_t m_trueTrackSize;
  std::size_t m_capturedTrackSize;

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

  /// True bound track parameters
  TVectorD* m_boundTrackParametersTruth = nullptr;
  TMatrixD* m_boundTrackCovTruth = nullptr;

  /// Initial guess of the momentum at the IP
  TLorentzVector* m_ipMomentumGuess = nullptr;

  /// Initial guess of the vertex at the IP
  TVector3* m_vertexGuess = nullptr;

  /// KF predicted momentum at the IP
  TLorentzVector* m_ipMomentumEst = nullptr;

  /// KF predicted vertex at the IP
  TVector3* m_vertexEst = nullptr;

  /// True momentum at the IP
  TLorentzVector* m_ipMomentumTruth = nullptr;

  /// True vertex at the IP
  TVector3* m_vertexTruth = nullptr;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
