#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"

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
#include "TrackingPipeline/Infrastructure/TypeDefinitions.hpp"
#include "TrackingPipeline/Io/E320RootDataReader.hpp"

/// @brief E320-specific track writer
class RootTrackWriter : public IWriter {
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
    /// Output tree name
    std::string treeName;
    /// Output file path
    std::string filePath;
  };

  RootTrackWriter(const RootTrackWriter &) = delete;
  RootTrackWriter(const RootTrackWriter &&) = delete;

  /// @brief constructor
  RootTrackWriter(const Config &config, Acts::Logging::Level level);

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

  std::unique_ptr<const Acts::Logger> m_logger;

  /// The output file
  TFile *m_file = nullptr;

  /// The output tree
  TTree *m_tree = nullptr;

 protected:
  /// Measurement hits in the surface frame
  std::vector<TVector2> m_trackHitsLocal;

  /// Measurement hits in the global frame
  std::vector<TVector3> m_trackHitsGlobal;

  /// Direction measurements in the track frame
  std::vector<TVector3> m_onSurfaceTrackDirection;

  /// Covariances of the track hits
  std::vector<TMatrixD> m_trackHitCovs;

  /// Covariances of the track agnles
  std::vector<TMatrixD> m_trackAngleCovs;

  /// Geometry ids of the track hits
  std::vector<std::size_t> m_geometryIds;

  /// KF predicted track hits in the surface frame
  std::vector<TVector2> m_predictedTrackHitsLocal;
  std::vector<TVector2> m_filteredTrackHitsLocal;
  std::vector<TVector2> m_smoothedTrackHitsLocal;

  /// KF predicted track hits in the global frame
  std::vector<TVector3> m_predictedTrackHitsGlobal;
  std::vector<TVector3> m_filteredTrackHitsGlobal;
  std::vector<TVector3> m_smoothedTrackHitsGlobal;

  /// KF predicted on surface momenta in the track frame
  std::vector<TLorentzVector> m_predictedOnSurfaceMomentum;
  std::vector<TLorentzVector> m_filteredOnSurfaceMomentum;
  std::vector<TLorentzVector> m_smoothedOnSurfaceMomentum;

  /// KF residuals with respect to the measurements
  std::vector<TVector2> m_predictedHitResiduals;
  std::vector<TVector2> m_filteredHitResiduals;
  std::vector<TVector2> m_smoothedHitResiduals;

  std::vector<TVector2> m_predictedAngleResiduals;
  std::vector<TVector2> m_filteredAngleResiduals;
  std::vector<TVector2> m_smoothedAngleResiduals;

  /// KF pulls with respect to the measurements
  std::vector<TVector2> m_predictedHitPulls;
  std::vector<TVector2> m_filteredHitPulls;
  std::vector<TVector2> m_smoothedHitPulls;

  std::vector<TVector2> m_predictedAnglePulls;
  std::vector<TVector2> m_filteredAnglePulls;
  std::vector<TVector2> m_smoothedAnglePulls;

  /// Guessed bound track parameters
  TVectorD m_boundTrackParametersGuess;
  TMatrixD m_boundTrackCovGuess;

  /// KF predicted bound track parameters
  TVectorD m_boundTrackParametersEst;
  TMatrixD m_boundTrackCovEst;

  /// Initial guess of the momentum at the IP
  TLorentzVector m_originMomentumGuess;

  /// Initial guess of the vertex at the IP
  TVector3 m_vertexGuess;

  /// KF predicted momentum at the IP
  TLorentzVector m_originMomentumEst;

  /// KF predicted vertex at the IP
  TVector3 m_vertexEst;

  /// Chi2 of the track
  /// with respect ot the
  /// measurement
  double m_chi2Predicted = 0;
  double m_chi2Filtered = 0;
  double m_chi2Smoothed = 0;

  /// Number of degrees of freedom
  /// of the track
  std::size_t m_ndf = 0;

  /// Track ID
  std::size_t m_trackId = 0;

  /// Event ID
  std::size_t m_eventId = 0;

  /// PDG ID
  int m_pdgId = 0;

  /// Charge
  int m_charge = 0;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
