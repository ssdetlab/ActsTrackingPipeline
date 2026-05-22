#pragma once

#include <memory>

#include "TFile.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TVector3.h"
#include "TrackingPipeline/Alignment/AlignmentStore.hpp"

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
  explicit AlignmentParametersProvider(const Config& config);

  /// Write out data to the input stream
  std::shared_ptr<AlignmentStore> getAlignmentStore();

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
  std::shared_ptr<AlignmentStore> m_store = nullptr;

  /// Detector element geometry ID
  std::vector<int>* m_geoId = nullptr;

  /// Detector element nominal transform
  std::vector<TVector3>* m_nominalTranslation = nullptr;
  std::vector<TMatrixD>* m_nominalRotation = nullptr;

  /// Detector element new aligned transform
  std::vector<TVector3>* m_newTranslation = nullptr;
  std::vector<TMatrixD>* m_newRotation = nullptr;

  std::size_t m_alignmentDof = 0;

  TMatrixD* m_alignmentCov = nullptr;
};
