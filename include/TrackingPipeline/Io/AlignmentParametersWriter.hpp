#pragma once

#include "Acts/Geometry/DetectorElementBase.hpp"

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
  using AlignmentParameters =
      std::unordered_map<Acts::DetectorElementBase *, Acts::Transform3>;

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

  ReadDataHandle<AlignmentParameters> m_alignmentResults{
      this, "InputAlignmentResults"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// The output file
  TFile *m_file = nullptr;

  /// The output tree
  TTree *m_tree = nullptr;

 protected:
  /// Detector element geometry ID
  int m_geoId;

  /// Detector element nominal transform
  TVector3 m_nominalTranslation;
  TMatrixD m_nominalRotation = TMatrixD(3, 3);

  /// Detector element new aligned transform
  TVector3 m_newTranslation;
  TMatrixD m_newRotation = TMatrixD(3, 3);

  /// Difference between nominal and aligned
  ///
  /// Note: to apply to the nominal the order is as follows:
  /// nominal.pretranslate(delta)
  /// nominal.rotate(delta)
  ///
  TVector3 m_deltaTranslation;
  TMatrixD m_deltaRotation = TMatrixD(3, 3);

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
