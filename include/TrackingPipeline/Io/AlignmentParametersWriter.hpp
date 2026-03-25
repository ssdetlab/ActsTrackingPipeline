#pragma once

#include "ActsAlignment/Kernel/Alignment.hpp"

#include "TFile.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TVector3.h"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

class AlignmentParametersWriter : public IWriter {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Input alignment results
    std::string inputAlignmentResults;
    /// Name of the input tree
    std::string treeName;
    /// The names of the input files
    std::string filePath;
  };

  AlignmentParametersWriter(const AlignmentParametersWriter &) = delete;
  AlignmentParametersWriter(const AlignmentParametersWriter &&) = delete;

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  AlignmentParametersWriter(const Config &config, Acts::Logging::Level level);

  /// @brief Finalize method
  ProcessCode finalize() override;

  /// Writer name() method
  std::string name() const override { return "AlignmentParametersWriter"; }

  /// Write out data to the input stream
  ProcessCode write(const AlgorithmContext &ctx) override;

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// Private access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  ReadDataHandle<ActsAlignment::AlignmentResult> m_alignmentResults{
      this, "InputAlignmentResults"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// The output file
  TFile *m_file = nullptr;

  /// The output tree
  TTree *m_tree = nullptr;

 protected:
  /// Detector element geometry ID
  std::vector<int> m_geoId;

  /// Detector element nominal transform
  std::vector<TVector3> m_nominalTranslation;
  std::vector<TMatrixD> m_nominalRotation;

  /// Detector element new aligned transform
  std::vector<TVector3> m_newTranslation;
  std::vector<TMatrixD> m_newRotation;

  std::size_t m_alignmentDof;

  TMatrixD m_alignmentCov;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
