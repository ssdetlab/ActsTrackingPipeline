#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/EventData/TrackParameters.hpp>

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
  struct Constraints {
    double minChi2;
    double maxChi2;

    double minVertexEstLong;
    double maxVertexEstLong;

    double minVertexEstShort;
    double maxVertexEstShort;

    double minAbsMomentumEst;
    double maxAbsMomentumEst;
  };

  /// @brief The nested configuration struct
  struct Config {
    /// Output source links
    std::string outputSourceLinks;
    /// Output seeds
    std::string outputSeedsGuess;
    /// Output track parameters guess
    std::string outputTrackParametersGuess;
    /// Output fitted seeds
    std::string outputSeedsEst;
    /// Output track parameters est
    std::string outputTrackParametersEst;
    /// The names of the input files
    std::vector<std::string> filePaths;
    /// Name of the input tree
    std::string treeName;
    /// Cuts
    Constraints constraints;
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
      this, "OutputSourceLinks"};

  /// WriteDataHandle for the seed data
  WriteDataHandle<IndexSeeds> m_outputSeedsGuess{this, "SeedsGuess"};

  WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_outputTrackParametersGuess{this, "OutputTrackParametersGuess"};

  /// WriteDataHandle for the fitted seed data
  WriteDataHandle<IndexSeeds> m_outputSeedsEst{this, "SeedsEst"};

  WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_outputTrackParametersEst{this, "OutputTrackParametersEst"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// Mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// Vector of {eventNr, entryMin, entryMax}
  std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> m_eventMap;

  /// The input tree name
  TFile* m_file = nullptr;
  TTree* m_tree = nullptr;
  TChain* m_chainOwner = nullptr;

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
  double m_chi2Predicted = 0;
  double m_chi2Filtered = 0;
  double m_chi2Smoothed = 0;

  /// Number of degrees of freedom
  /// of the track
  std::size_t m_ndf = 0;

  /// TrackId
  std::size_t m_trackId = 0;

  /// EventId
  std::size_t m_eventId = 0;

  /// PDG ID
  int m_pdgId = 0;

  /// Charge
  int m_charge = 0;

  /// Guessed bound track parameters
  TVectorD* m_boundTrackParametersGuess = nullptr;
  TMatrixD* m_boundTrackCovGuess = nullptr;

  /// KF predicted bound track parameters
  TVectorD* m_boundTrackParametersEst = nullptr;
  TMatrixD* m_boundTrackCovEst = nullptr;

  /// Initial guess of the momentum at the IP
  TLorentzVector* m_originMomentumGuess = nullptr;

  /// Initial guess of the vertex at the IP
  TVector3* m_vertexGuess = nullptr;

  /// KF predicted momentum at the IP
  TLorentzVector* m_originMomentumEst = nullptr;

  /// KF predicted vertex at the IP
  TVector3* m_vertexEst = nullptr;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
