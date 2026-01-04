#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TVector3.h"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

using namespace Acts::UnitLiterals;

using TrackID = std::tuple<int, int, int>;

/// @brief Writer to store fitted track data in
/// ROOT file
///
/// Writer that accepts fitted track data from KF
/// derives the basic performance metrics, such as
/// chi2 and residuals, and stores them in a ROOT file.
///
/// @note Assumes that the tracks are simulated and
/// the truth information is available
class RootSimTrackCandidateWriter : public IWriter {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Surface accessor
    Acts::SourceLinkSurfaceAccessor surfaceAccessor;
    /// Fitted track collection
    std::string inputTrackCandidates;
    /// Truth cluster data
    std::string inputTruthClusters;
    /// Name of the input tree
    std::string treeName;
    /// The names of the input files
    std::string filePath;
  };

  RootSimTrackCandidateWriter(const RootSimTrackCandidateWriter &) = delete;
  RootSimTrackCandidateWriter(const RootSimTrackCandidateWriter &&) = delete;

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  RootSimTrackCandidateWriter(const Config &config, Acts::Logging::Level level);

  /// @brief Finalize method
  ProcessCode finalize() override;

  /// Writer name() method
  std::string name() const override { return "RootSimTrackCandidateWriter"; }

  /// Write out data to the input stream
  ProcessCode write(const AlgorithmContext &ctx) override;

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// Private access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  ReadDataHandle<Tracks> m_inputTrackCandidates{this, "InputTrackCandidates"};

  ReadDataHandle<SimClusters> m_inputTruthClusters{this, "InputTruthClusters"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// The output file
  TFile *m_file = nullptr;

  /// The output tree
  TTree *m_tree = nullptr;

 protected:
  /// True hits
  std::vector<TVector3> m_trueTrackHitsGlobal;
  std::vector<TVector2> m_trueTrackHitsLocal;
  std::vector<TLorentzVector> m_onSurfaceMomentum;
  std::vector<int> m_isSignal;

  /// Measurement hits
  std::vector<TVector3> m_trackHitsGlobal;
  std::vector<TVector2> m_trackHitsLocal;

  /// Covariances of the track hits
  std::vector<TMatrixD> m_trackHitCovs;

  /// Geometry ids of the track hits
  std::vector<std::size_t> m_geometryIds;

  /// KF predicted track hits
  std::vector<TVector3> m_predictedTrackHitsGlobal;
  std::vector<TVector3> m_filteredTrackHitsGlobal;

  std::vector<TVector2> m_predictedTrackHitsLocal;
  std::vector<TVector2> m_filteredTrackHitsLocal;

  /// KF residuals with respect to the true hits
  std::vector<TVector2> m_truePredictedResiduals;
  std::vector<TVector2> m_trueFilteredResiduals;

  /// KF residuals with respect to the measurements
  std::vector<TVector2> m_predictedResiduals;
  std::vector<TVector2> m_filteredResiduals;

  /// KF pulls with respect to the true hits
  std::vector<TVector2> m_truePredictedPulls;
  std::vector<TVector2> m_trueFilteredPulls;

  /// KF pulls with respect to the measurements
  std::vector<TVector2> m_predictedPulls;
  std::vector<TVector2> m_filteredPulls;

  /// Chi2 of the track
  /// with respect ot the
  /// measurement
  double m_chi2Predicted;
  double m_chi2Filtered;

  /// Number of degrees of freedom
  /// of the track
  std::size_t m_ndf;

  /// TrackId
  std::vector<std::size_t> m_stateTrackId;
  std::vector<std::size_t> m_stateParentTrackId;
  std::vector<std::size_t> m_stateRunId;

  std::size_t m_trackId;
  std::size_t m_parentTrackId;
  std::size_t m_runId;

  /// EventId
  std::size_t m_eventId;

  /// True track size
  std::size_t m_trueTrackSize;
  std::size_t m_trackInCandidateSize;

  /// PDG ID
  int m_pdgId;

  /// Charge
  int m_charge;

  /// Initial guess of the momentum at the IP
  TLorentzVector m_ipMomentumGuess;
  TVector3 m_vertexGuess;

  /// True momentum at the IP
  TLorentzVector m_ipMomentumTruth;
  TVector3 m_vertexTruth;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
