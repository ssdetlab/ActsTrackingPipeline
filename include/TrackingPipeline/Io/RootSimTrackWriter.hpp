#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/EventData/TrackParameters.hpp>

#include <cstddef>
#include <vector>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

using namespace Acts::UnitLiterals;

using TrackID = std::tuple<int, int, int>;

/// @brief writer storing fitted sim tracks in a ROOT file
class RootSimTrackWriter : public IWriter {
 public:
  /// @brief nested configuration struct
  struct Config {
    /// Surface accessor
    Acts::SourceLinkSurfaceAccessor surfaceAccessor;
    /// Reference surface
    const Acts::Surface *referenceSurface;
    /// Input track container
    std::string inputTrackContainer;
    /// Input track index containers
    std::string inputTracks;
    /// Input track container
    std::string inputTrackParametersGuesses;
    /// Input sim clusters
    std::string inputSimClusters;
    /// Output tree name
    std::string treeName;
    /// Output file path
    std::string filePath;
  };

  RootSimTrackWriter(const RootSimTrackWriter &) = delete;
  RootSimTrackWriter(const RootSimTrackWriter &&) = delete;

  /// @brief constructor
  RootSimTrackWriter(const Config &config, Acts::Logging::Level level);

  /// @brief finalize method
  ProcessCode finalize() override;

  /// @brief get writer name
  std::string name() const override { return "RootFittedTrackWriter"; }

  /// @brief write out data to the file
  ProcessCode write(const AlgorithmContext &ctx) override;

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// Private access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  ReadDataHandle<KFFitterTrackContainer> m_inputTrackContainer{
      this, "InputTrackContainer"};

  ReadDataHandle<IndexTracks> m_inputTracks{this, "InputTracks"};

  ReadDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_inputTrackParametersGuesses{this, "InputTrackParametersGuesses"};

  ReadDataHandle<SimClusters> m_inputSimClusters{this, "InputSimClusters"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// The output file
  TFile *m_file = nullptr;

  /// The output tree
  TTree *m_tree = nullptr;

 protected:
  /// True hits
  std::vector<TVector3> m_trueTrackHitsGlobal;
  std::vector<TVector2> m_trueTrackHitsLocal;
  std::vector<TLorentzVector> m_onSurfaceMomentumTruth;
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
  std::vector<TVector3> m_smoothedTrackHitsGlobal;

  std::vector<TVector2> m_predictedTrackHitsLocal;
  std::vector<TVector2> m_filteredTrackHitsLocal;
  std::vector<TVector2> m_smoothedTrackHitsLocal;

  /// KF residuals with respect to the true hits
  std::vector<TVector2> m_truePredictedResiduals;
  std::vector<TVector2> m_trueFilteredResiduals;
  std::vector<TVector2> m_trueSmoothedResiduals;

  /// KF residuals with respect to the measurements
  std::vector<TVector2> m_predictedResiduals;
  std::vector<TVector2> m_filteredResiduals;
  std::vector<TVector2> m_smoothedResiduals;

  /// KF pulls with respect to the true hits
  std::vector<TVector2> m_truePredictedPulls;
  std::vector<TVector2> m_trueFilteredPulls;
  std::vector<TVector2> m_trueSmoothedPulls;

  /// KF pulls with respect to the measurements
  std::vector<TVector2> m_predictedPulls;
  std::vector<TVector2> m_filteredPulls;
  std::vector<TVector2> m_smoothedPulls;

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
  std::vector<std::size_t> m_stateTrackId;
  std::vector<std::size_t> m_stateParentTrackId;
  std::vector<std::size_t> m_stateRunId;

  std::size_t m_trackId = 0;
  std::size_t m_parentTrackId = 0;
  std::size_t m_runId = 0;

  /// EventId
  std::size_t m_eventId = 0;

  /// True track size
  std::size_t m_trueTrackSize = 0;
  std::size_t m_capturedTrackSize = 0;

  /// PDG ID
  int m_pdgId = 0;

  /// Charge
  int m_charge = 0;

  /// Guessed bound track parameters
  TVectorD m_boundTrackParametersGuess;
  TMatrixD m_boundTrackCovGuess;

  /// KF predicted bound track parameters
  TVectorD m_boundTrackParametersEst;
  TMatrixD m_boundTrackCovEst;

  /// True bound track parameters
  TVectorD m_boundTrackParametersTruth;
  TMatrixD m_boundTrackCovTruth;

  /// Initial guess of the momentum at the IP
  TLorentzVector m_originMomentumGuess;

  /// Initial guess of the vertex at the IP
  TVector3 m_vertexGuess;

  /// KF predicted momentum at the IP
  TLorentzVector m_originMomentumEst;

  /// KF predicted vertex at the IP
  TVector3 m_vertexEst;

  /// True momentum at the IP
  TLorentzVector m_originMomentumTruth;

  /// True vertex at the IP
  TVector3 m_vertexTruth;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
