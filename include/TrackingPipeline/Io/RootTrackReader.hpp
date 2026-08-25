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
#include "TrackingPipeline/EventData/ExtendedSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

class RootTrackReader : public IReader {
 public:
  /// Exntended sizes shorthands
  static constexpr const std::size_t ExtendedLocalSize =
      ExtendedSourceLink::localSubspaceSize;
  static constexpr const std::size_t ExtendedGlobalSize =
      ExtendedSourceLink::globalSubspaceSize;

  /// @brief Cuts on the tracks read
  struct Constraints {
    /// Track smoothed chi2 range
    double minSmoothedChi2;
    double maxSmoothedChi2;
    /// Transverse vertex range
    double minVertexEstLong;
    double maxVertexEstLong;
    double minVertexEstShort;
    double maxVertexEstShort;
    /// Momentum range
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
    /// Cuts on the read tracks
    Constraints constraints;
    /// Merge into a single event flag
    bool mergeIntoOneEvent;
    /// Backwards propagation flag
    bool backwards;
  };

  RootTrackReader(const RootTrackReader&) = delete;
  RootTrackReader(const RootTrackReader&&) = delete;

  /// @brief Constructor
  ///
  /// @param config configuration struct
  /// @param level logging level
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

  /// WriteDataHandle for the source links data
  WriteDataHandle<std::vector<Acts::SourceLink>> m_outputSourceLinks{
      this, "OutputSourceLinks"};

  /// WriteDataHandle for the guess seed data
  WriteDataHandle<IndexSeeds> m_outputSeedsGuess{this, "SeedsGuess"};

  WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_outputTrackParametersGuess{this, "OutputTrackParametersGuess"};

  /// WriteDataHandle for the fitted seed data
  WriteDataHandle<IndexSeeds> m_outputSeedsEst{this, "SeedsEst"};

  WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_outputTrackParametersEst{this, "OutputTrackParametersEst"};

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// Vector of {eventNr, entryMin, entryMax}
  std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> m_eventMap;

  /// The IO handles
  TFile* m_file = nullptr;
  TTree* m_tree = nullptr;
  TChain* m_chainOwner = nullptr;

  /// List of constraint surface IDs
  std::unordered_set<std::size_t> m_constraintSurfacesGeoIds;

 protected:
  /// Measurement hits in the surface frame
  std::vector<TVector2>* m_trackHitsLocal = nullptr;

  /// Measurement hits in the global frame
  std::vector<TVector3>* m_trackHitsGlobal = nullptr;

  /// Direction measurements in the track frame
  std::vector<TVector3>* m_onSurfaceTrackDirection = nullptr;

  /// Covariances of the track hits
  std::vector<TMatrixD>* m_trackHitCovs = nullptr;

  /// Covariances of the track agnles
  std::vector<TMatrixD>* m_trackAngleCovs = nullptr;

  /// Geometry ids of the track hits
  std::vector<std::size_t>* m_geometryIds = nullptr;

  /// KF predicted track hits in the surface frame
  std::vector<TVector2>* m_predictedTrackHitsLocal = nullptr;
  std::vector<TVector2>* m_filteredTrackHitsLocal = nullptr;
  std::vector<TVector2>* m_smoothedTrackHitsLocal = nullptr;

  /// KF predicted track hits in the global frame
  std::vector<TVector3>* m_predictedTrackHitsGlobal = nullptr;
  std::vector<TVector3>* m_filteredTrackHitsGlobal = nullptr;
  std::vector<TVector3>* m_smoothedTrackHitsGlobal = nullptr;

  /// KF predicted on surface momenta in the track frame
  std::vector<TLorentzVector>* m_predictedOnSurfaceMomentum = nullptr;
  std::vector<TLorentzVector>* m_filteredOnSurfaceMomentum = nullptr;
  std::vector<TLorentzVector>* m_smoothedOnSurfaceMomentum = nullptr;

  /// KF residuals with respect to the measurements
  std::vector<TVector2>* m_predictedHitResiduals = nullptr;
  std::vector<TVector2>* m_filteredHitResiduals = nullptr;
  std::vector<TVector2>* m_smoothedHitResiduals = nullptr;

  std::vector<TVector2>* m_predictedAngleResiduals = nullptr;
  std::vector<TVector2>* m_filteredAngleResiduals = nullptr;
  std::vector<TVector2>* m_smoothedAngleResiduals = nullptr;

  /// KF pulls with respect to the measurements
  std::vector<TVector2>* m_predictedHitPulls = nullptr;
  std::vector<TVector2>* m_filteredHitPulls = nullptr;
  std::vector<TVector2>* m_smoothedHitPulls = nullptr;

  std::vector<TVector2>* m_predictedAnglePulls = nullptr;
  std::vector<TVector2>* m_filteredAnglePulls = nullptr;
  std::vector<TVector2>* m_smoothedAnglePulls = nullptr;

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
