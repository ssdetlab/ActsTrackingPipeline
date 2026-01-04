#pragma once

#include <memory>

#include "TFile.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TVector3.h"
#include "TrackingPipeline/Alignment/AlignmentContext.hpp"

class AlignmentParametersProvider {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Name of the input tree
    std::string treeName;
    /// The names of the input files
    std::string filePath;
  };

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  AlignmentParametersProvider(const Config& config);

  /// Write out data to the input stream
  std::shared_ptr<AlignmentContext::AlignmentStore> getAlignmentStore();

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// The config class
  Config m_cfg;

  /// The input file
  TFile* m_file = nullptr;

  /// The input tree
  TTree* m_tree = nullptr;

 protected:
  /// Alignment store
  std::shared_ptr<AlignmentContext::AlignmentStore> m_store = nullptr;

  /// Detector element geometry ID
  int m_geoId;

  /// Detector element nominal transform
  TVector3* m_nominalTranslation = nullptr;
  TMatrixD* m_nominalRotation = nullptr;

  /// Detector element new aligned transform
  TVector3* m_newTranslation = nullptr;
  TMatrixD* m_newRotation = nullptr;

  /// Difference between nominal and aligned
  ///
  /// Note: to apply to the nominal the order is as follows:
  /// nominal.pretranslate(delta)
  /// nominal.rotate(delta)
  ///
  TVector3* m_deltaTranslation = nullptr;
  TMatrixD* m_deltaRotation = nullptr;
};
