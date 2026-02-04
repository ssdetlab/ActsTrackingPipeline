#pragma once

#include "Acts/Utilities/Logger.hpp"

#include "TFile.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVectorD.h"
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
class RootSimClusterWriter : public IWriter {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Truth cluster data
    std::string inputClusters;
    /// Name of the input tree
    std::string treeName;
    /// The names of the input files
    std::string filePath;
  };

  RootSimClusterWriter(const RootSimClusterWriter &) = delete;
  RootSimClusterWriter(const RootSimClusterWriter &&) = delete;

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  RootSimClusterWriter(const Config &config, Acts::Logging::Level level);

  /// @brief Finalize method
  ProcessCode finalize() override;

  /// Writer name() method
  std::string name() const override { return "RootSimClusterWriter"; }

  /// Write out data to the input stream
  ProcessCode write(const AlgorithmContext &ctx) override;

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// Private access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  ReadDataHandle<SimClusters> m_inputClusters{this, "InputClusters"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// The output file
  TFile *m_file = nullptr;

  /// The output tree
  TTree *m_tree = nullptr;

 protected:
  /// Cluster parameters
  TVector3 m_geoCenterGlobal;
  TVector2 m_geoCenterLocal;
  TMatrixD m_clusterCov = TMatrixD(2, 2);

  std::size_t m_geoId;
  std::size_t m_eventId;

  int m_isSignal;

  /// Measurement hits
  std::vector<TVector3> m_trackHitsGlobal;
  std::vector<TVector2> m_trackHitsLocal;

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
  std::vector<TLorentzVector> m_onSurfaceMomentum;

  std::vector<int> m_charge;
  std::vector<int> m_pdgId;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
