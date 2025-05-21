#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

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

/// @brief Writer to store fitted track data in
/// ROOT file
///
/// Writer that accepts fitted track data from KF
/// derives the basic performance metrics, such as
/// chi2 and residuals, and stores them in a ROOT file.
class RootTrackWriter : public IWriter {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Surface accessor
    Acts::SourceLinkSurfaceAccessor surfaceAccessor;
    /// Fitted track collection
    std::string inputTracks;
    /// Name of the input tree
    std::string treeName;
    /// The names of the input files
    std::string filePath;
  };

  RootTrackWriter(const RootTrackWriter &) = delete;
  RootTrackWriter(const RootTrackWriter &&) = delete;

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  RootTrackWriter(const Config &config, Acts::Logging::Level level);

  /// @brief Finalize method
  ProcessCode finalize() override;

  /// Writer name() method
  std::string name() const override { return "RootTrackWriter"; }

  /// Write out data to the input stream
  ProcessCode write(const AlgorithmContext &ctx) override;

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// Private access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  ReadDataHandle<Tracks> m_inputTracks{this, "Tracks"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// The output file
  TFile *m_file = nullptr;

  /// The output tree
  TTree *m_tree = nullptr;

 protected:
  /// Measurement hits
  std::vector<TVector3> m_trackHitsGlobal;
  std::vector<TVector2> m_trackHitsLocal;

  /// Covariances of the track hits
  std::vector<TMatrixD> m_trackHitCovs;

  /// Geometry ids of the track hits
  std::vector<int> m_geometryIds;

  /// KF predicted track hits
  std::vector<TVector3> m_predictedTrackHitsGlobal;
  std::vector<TVector3> m_filteredTrackHitsGlobal;
  std::vector<TVector3> m_smoothedTrackHitsGlobal;

  std::vector<TVector2> m_predictedTrackHitsLocal;
  std::vector<TVector2> m_filteredTrackHitsLocal;
  std::vector<TVector2> m_smoothedTrackHitsLocal;

  /// KF residuals with respect to the measurements
  std::vector<TVector2> m_predictedResiduals;
  std::vector<TVector2> m_filteredResiduals;
  std::vector<TVector2> m_smoothedResiduals;

  /// KF pulls with respect to the measurements
  std::vector<TVector2> m_predictedPulls;
  std::vector<TVector2> m_filteredPulls;
  std::vector<TVector2> m_smoothedPulls;

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
  std::vector<int> m_trackId;

  /// EventId
  int m_eventId;

  /// Initial guess of the momentum at the IP
  TLorentzVector m_ipMomentumGuess;
  TVector3 m_vertexGuess;

  /// KF predicted momentum at the IP
  TLorentzVector m_ipMomentumEst;
  TVector3 m_ipMomentumError;
  TVector3 m_vertexEst;
  TVector3 m_vertexError;

  /// PDG ID
  int m_pdgId;

  /// Charge
  int m_charge;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
