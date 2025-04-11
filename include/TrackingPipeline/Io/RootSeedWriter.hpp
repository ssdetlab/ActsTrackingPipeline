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

class RootSeedWriter : public IWriter {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Fitted track collection
    std::string inputSeeds;
    /// Name of the input tree
    std::string treeName;
    /// The names of the input files
    std::string filePath;
  };

  RootSeedWriter(const RootSeedWriter &) = delete;
  RootSeedWriter(const RootSeedWriter &&) = delete;

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  RootSeedWriter(const Config &config, Acts::Logging::Level level);

  /// @brief Finalize method
  ProcessCode finalize() override;

  /// Writer name() method
  std::string name() const override { return "RootSeedWriter"; }

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

  std::unique_ptr<const Acts::Logger> m_logger;

  /// The output file
  TFile *m_file = nullptr;

  /// The output tree
  TTree *m_tree = nullptr;

 protected:
  /// Number of source links
  /// in a seed
  int m_size;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
