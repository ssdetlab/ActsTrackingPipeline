#pragma once

#include <cstddef>

#include <TLorentzVector.h>

#include "TFile.h"
#include "TTree.h"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

using TrackID = std::tuple<std::int32_t, std::int32_t, std::int32_t>;

class RootSimSeedWriter : public IWriter {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Fitted track collection
    std::string inputSeeds;
    /// Truth cluster data
    std::string inputTruthClusters;
    /// Name of the input tree
    std::string treeName;
    /// The names of the input files
    std::string filePath;
    /// Target size of the true track
    std::size_t targetTrueTrackSize;
  };

  RootSimSeedWriter(const RootSimSeedWriter &) = delete;
  RootSimSeedWriter(const RootSimSeedWriter &&) = delete;

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  RootSimSeedWriter(const Config &config, Acts::Logging::Level level);

  /// @brief Finalize method
  ProcessCode finalize() override;

  /// Writer name() method
  std::string name() const override { return "RootSimSeedWriter"; }

  /// Write out data to the input stream
  ProcessCode write(const AlgorithmContext &ctx) override;

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// Private access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  ReadDataHandle<Seeds> m_seeds{this, "InputSeeds"};

  ReadDataHandle<SimClusters> m_truthClusters{this, "TruthClusters"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// The output file
  TFile *m_file = nullptr;

  /// The output tree
  TTree *m_tree = nullptr;

 protected:
  /// Fraction of a true track
  /// contained within the seed
  /// source link list
  double m_efficiency;

  /// Number of source links
  /// in a seed
  int m_size;

  /// Flag idicating if pivot
  /// cluster is a signal
  bool m_isSignal;

  /// True momentum at the IP
  /// associated with a pivot
  TLorentzVector m_ipMomentumTruth;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
