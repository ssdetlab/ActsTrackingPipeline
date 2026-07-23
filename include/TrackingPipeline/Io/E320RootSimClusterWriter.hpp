#pragma once

#include "Acts/Utilities/Logger.hpp"
#include <Acts/EventData/SourceLink.hpp>

#include <cstddef>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/ExtendedSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

namespace E320 {

using namespace Acts::UnitLiterals;

using TrackID = std::tuple<int, int, int>;

/// @brief RootSimClusterWriter specialization implementing
/// Detector-upstream constraints split
class E320RootSimClusterWriter : public IWriter {
 public:
  /// Exntended sizes shorthands
  static constexpr const std::size_t ExtendedLocalSize =
      ExtendedSourceLink::localSubspaceSize;
  static constexpr const std::size_t ExtendedGlobalSize =
      ExtendedSourceLink::globalSubspaceSize;

  /// @brief The nested configuration struct
  struct Config {
    /// Input sim clusters
    std::string inputClusters;
    /// Output tree name
    std::string treeName;
    /// Output file path
    std::string filePath;
  };

  E320RootSimClusterWriter(const E320RootSimClusterWriter &) = delete;
  E320RootSimClusterWriter(const E320RootSimClusterWriter &&) = delete;

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  E320RootSimClusterWriter(const Config &config, Acts::Logging::Level level);

  /// @brief Finalize method
  ProcessCode finalize() override;

  /// Writer name() method
  std::string name() const override { return "E320RootSimClusterWriter"; }

  /// Write out data to the input stream
  ProcessCode write(const AlgorithmContext &ctx) override;

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// Private access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  /// ReadDataHandles of the data
  ReadDataHandle<SimClusters> m_inputSimClusters{this, "SimClusters"};

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// The output file
  TFile *m_file = nullptr;

  /// The output tree
  TTree *m_tree = nullptr;

 protected:
  /// Cluster hit in the surface frame
  TVector2 m_geoCenterLocal;

  /// Cluster hit in the global frame
  TVector3 m_geoCenterGlobal;

  /// Cluster hit covariance in the surface frame
  TMatrixD m_clusterCov = TMatrixD(2, 2);

  /// Momentum direction measurement in the track frame
  TVector3 m_onSurfaceDirection;

  /// Angle covariance in the surface frame
  TMatrixD m_angleCov = TMatrixD(2, 2);

  /// Surface geometry ID
  std::size_t m_geoId = 0;

  /// Event ID
  std::size_t m_eventId = 0;

  /// Signal/background flag
  int m_isSignal = 0;

  /// True track hits within the clusters
  std::vector<TVector3> m_trackHitsGlobal;
  std::vector<TVector2> m_trackHitsLocal;

  /// Track IDs
  std::vector<int> m_trackId;
  std::vector<int> m_parentTrackId;
  std::vector<int> m_runId;

  /// Bound origin parameters
  std::vector<TVectorD> m_boundTrackParameters;
  std::vector<TMatrixD> m_boundTrackCov;

  /// Origin momentum
  std::vector<TLorentzVector> m_originMomentum;

  /// Origin vertex
  std::vector<TVector3> m_vertex;

  /// Momentum at clusters
  std::vector<TLorentzVector> m_onSurfaceMomentumTruth;

  /// Charges of the tracks' particles
  std::vector<int> m_charge;

  /// PDG IDs of the tracks' particles
  std::vector<int> m_pdgId;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};

}  // namespace E320
