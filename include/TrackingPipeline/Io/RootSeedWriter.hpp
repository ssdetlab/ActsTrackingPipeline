
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

/// @brief writer storing sim seeds in a ROOT file
class RootSeedWriter : public IWriter {
 public:
  /// @brief nested configuration struct
  struct Config {
    /// Input seeds
    std::string inputSeeds;
    /// Input track parameters
    std::string inputTrackParameters;
    /// Input source link
    std::string inputSourceLinks;
    /// Outout tree name
    std::string treeName;
    /// Output file name
    std::string filePath;
  };

  RootSeedWriter(const RootSeedWriter &) = delete;
  RootSeedWriter(const RootSeedWriter &&) = delete;

  /// @brief constructor
  RootSeedWriter(const Config &config, Acts::Logging::Level level);

  /// @brief finalize method
  ProcessCode finalize() override;

  /// @brief get the writer name
  std::string name() const override { return "RootSeedWriter"; }

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

  /// Estimated momentum, vertex at the IP
  /// associated with a pivot
  TLorentzVector m_originMomentumEst;
  TVector3 m_vertexEst;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
