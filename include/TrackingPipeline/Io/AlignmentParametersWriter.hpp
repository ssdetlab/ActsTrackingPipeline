#pragma once

#include "TFile.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TVector3.h"
#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

/// @brief writer of the alignment parameters from the alignment algorithm
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
  /// @param level logging level
  AlignmentParametersWriter(const Config &config, Acts::Logging::Level level);

  /// @brief Finalize method
  ProcessCode finalize() override;

  /// @brief Writer name() method
  std::string name() const override { return "AlignmentParametersWriter"; }

  /// @brief Write out data to the input stream
  ///
  /// @param ctx current algorithm context
  ///
  /// @return algorithm process code
  ProcessCode write(const AlgorithmContext &ctx) override;

  /// @brief Readonly access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// @brief Private access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  /// Configuration
  Config m_cfg;

  /// Read data handle
  ReadDataHandle<AlignmentAlgorithm::AlignmentResult> m_alignmentResults{
      this, "InputAlignmentResults"};

  /// Logging instance
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

  /// Number of alignment d.o.fs
  std::size_t m_alignmentDof = 0;

  /// Alignment parameters covariance
  TMatrixD m_alignmentCov;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
