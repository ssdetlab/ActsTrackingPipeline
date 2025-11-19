#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>

#include <TMatrixDfwd.h>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TVector3.h"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

class RootMeasurementWriter : public IWriter {
 public:
  using HitData = std::tuple<TVector3, TVector3, TLorentzVector>;

  /// @brief The nested configuration struct
  struct Config {
    /// Measurement data
    std::string inputMeasurements;
    /// Name of the output tree
    std::string treeName;
    /// The names of the output file
    std::string filePath;
  };

  RootMeasurementWriter(const RootMeasurementWriter &) = delete;
  RootMeasurementWriter(const RootMeasurementWriter &&) = delete;

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  RootMeasurementWriter(const Config &config, Acts::Logging::Level level);

  /// @brief Finalize method
  ProcessCode finalize() override;

  /// Writer name() method
  std::string name() const override { return "RootMeasurementWriter"; }

  /// Write out data to the input stream
  ProcessCode write(const AlgorithmContext &ctx) override;

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// Private access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  ReadDataHandle<std::vector<Acts::SourceLink>> m_inputMeasurements{
      this, "inputMeasurements"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// The output file
  TFile *m_file = nullptr;

  /// The output tree
  TTree *m_tree = nullptr;

 protected:
  TVector3 m_geoCenterGlobal;
  TVector2 m_geoCenterLocal;
  TMatrixD m_cov = TMatrixD(2, 2);

  std::size_t m_geoId;
  std::size_t m_eventId;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
