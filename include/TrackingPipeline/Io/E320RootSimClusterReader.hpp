#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <unordered_map>

#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVector3.h"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/ExtendedSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

namespace E320 {

/// @brief RootSimClusterReader specialization implementing
/// Detector-upstream constraints split
class E320RootSimClusterReader : public IReader {
 public:
  /// Exntended sizes shorthands
  static constexpr const std::size_t ExtendedLocalSize =
      ExtendedSourceLink::localSubspaceSize;
  static constexpr const std::size_t ExtendedGlobalSize =
      ExtendedSourceLink::globalSubspaceSize;

  /// @brief The nested configuration struct
  struct Config {
    /// Output source links
    std::string outputSourceLinks;
    /// Output sim clusters
    std::string outputSimClusters;
    /// Output detector source link indices
    std::string outputDetSourceLinkIndices;
    /// Output BPM source link indices
    std::string outputConstraintSourceLinkIndices;
    /// The names of the input files
    std::vector<std::string> filePaths;
    /// Name of the input tree
    std::string treeName;
    /// Geometry ID scope
    int minGeoId;
    int maxGeoId;
    /// Wheter to employ surfaces for local to global conversion
    bool surfaceLocalToGlobal;
    /// Surface map for high-precision local to global conversion
    std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
        surfaceMap;
    /// Backwards propagation flag
    bool backwards;
  };

  E320RootSimClusterReader(const E320RootSimClusterReader&) = delete;
  E320RootSimClusterReader(const E320RootSimClusterReader&&) = delete;

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  /// @param level Logging level
  E320RootSimClusterReader(const Config& config, Acts::Logging::Level level);

  /// Writer name() method
  std::string name() const override { return "E320RootSimClusterReader"; }

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

  /// WriteDataHandles of the data
  WriteDataHandle<std::vector<Acts::SourceLink>> m_outputSourceLinks{
      this, "OutputSourceLinks"};

  WriteDataHandle<SimClusters> m_outputSimClusters{this, "OutputSimClusters"};

  WriteDataHandle<std::vector<std::size_t>> m_outputDetSourceLinkIndices{
      this, "OutputDetSourceLinkIndices"};

  WriteDataHandle<std::vector<std::size_t>> m_outputConstraintSourceLinkIndices{
      this, "OutputConstraintSourceLinkIndices"};

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// Vector of {eventNr, entryMin, entryMax}
  std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> m_eventMap;

  /// The input tree name
  TTree* m_tree = nullptr;
  TFile* m_file = nullptr;
  TChain* m_chainOwner = nullptr;

  /// List of constraint surface IDs
  std::unordered_set<std::size_t> m_constraintSurfacesGeoIds;

 protected:
  /// Cluster hit in the surface frame
  TVector2* m_geoCenterLocal = nullptr;

  /// Cluster hit in the global frame
  TVector3* m_geoCenterGlobal = nullptr;

  /// Cluster hit covariance in the surface frame
  TMatrixD* m_clusterCov = nullptr;

  /// Momentum direction measurement in the track frame
  TVector3* m_onSurfaceDirection = nullptr;

  /// Angle covariance in the surface frame
  TMatrixD* m_angleCov = nullptr;

  /// Surface geometry ID
  std::size_t m_geoId = 0;

  /// Event ID
  std::size_t m_eventId = 0;

  /// Signal/background flag
  int m_isSignal = 0;

  /// True track hits within the clusters
  std::vector<TVector2>* m_trackHitsLocal = nullptr;
  std::vector<TVector3>* m_trackHitsGlobal = nullptr;

  /// Track IDs
  std::vector<int>* m_trackId = nullptr;
  std::vector<int>* m_parentTrackId = nullptr;
  std::vector<int>* m_runId = nullptr;

  /// Bound origin parameters
  std::vector<TVectorD>* m_boundTrackParameters = nullptr;
  std::vector<TMatrixD>* m_boundTrackCov = nullptr;

  /// Origin momentum
  std::vector<TLorentzVector>* m_originMomentum = nullptr;

  /// Origin vertex
  std::vector<TVector3>* m_vertex = nullptr;

  /// Momentum at clusters
  std::vector<TLorentzVector>* m_onSurfaceMomentumTruth = nullptr;

  /// Charges of the tracks' particles
  std::vector<int>* m_charge = nullptr;

  /// PDG IDs of the tracks' particles
  std::vector<int>* m_pdgId = nullptr;
};

}  // namespace E320
