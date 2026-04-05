#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include <Acts/EventData/TrackParameters.hpp>

#include <cstddef>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

using TrackID = std::tuple<int, int, int>;

/// @brief writer storing sim seeds in a ROOT file
class RootSimSeedWriter : public IWriter {
 public:
  /// @brief nested configuration struct
  struct Config {
    /// Input seeds
    std::string inputSeeds;
    /// Input track parameters
    std::string inputTrackParameters;
    /// Input source link
    std::string inputSourceLinks;
    /// Input sim cluster
    std::string inputSimClusters;
    /// Outout tree name
    std::string treeName;
    /// Output file name
    std::string filePath;
  };

  RootSimSeedWriter(const RootSimSeedWriter &) = delete;
  RootSimSeedWriter(const RootSimSeedWriter &&) = delete;

  /// @brief constructor
  RootSimSeedWriter(const Config &config, Acts::Logging::Level level);

  /// @brief finalize method
  ProcessCode finalize() override;

  /// @brief get the writer name
  std::string name() const override { return "RootSimSeedWriter"; }

  /// @brief write out data to the file
  ProcessCode write(const AlgorithmContext &ctx) override;

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

 protected:
  /// @brief access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

 private:
  /// The config class
  Config m_cfg;

  ReadDataHandle<IndexSeeds> m_inputSeeds{this, "InputSeeds"};

  ReadDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_inputTrackParameters{this, "InputTrackParameters"};

  ReadDataHandle<SimClusters> m_inputSimClusters{this, "InputSimClusters"};

  ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
      this, "InputSourceLinks"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// The output file
  TFile *m_file = nullptr;

  /// The output tree
  TTree *m_tree = nullptr;

 protected:
  /// Measurements in a given seed
  std::vector<TVector3> m_seedMeasurementsGlob;
  std::vector<TVector2> m_seedMeasurementsLoc;
  std::vector<int> m_geoIds;

  /// Event ID of the seed
  std::size_t m_eventId = 0;

  /// Number of source links
  /// in a seed
  std::size_t m_size = 0;

  /// Track ID of the seed
  std::size_t m_trackId = 0;

  /// Parent track ID of the seed
  std::size_t m_parentTrackId = 0;

  /// Run ID of the seed
  std::size_t m_runId = 0;

  /// Size of the true track
  std::size_t m_trueTrackSize = 0;

  /// Size of the true track in seed
  std::size_t m_trackInSeedSize = 0;

  /// Flag idicating if pivot
  /// cluster is a signal
  bool m_isSignal = false;

  /// True momentum, vertex at the IP
  /// associated with a pivot
  TLorentzVector m_originMomentumTruth;
  TVector3 m_vertexTruth;

  /// Estimated momentum, vertex at the IP
  /// associated with a pivot
  TLorentzVector m_originMomentumEst;
  TVector3 m_vertexEst;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
