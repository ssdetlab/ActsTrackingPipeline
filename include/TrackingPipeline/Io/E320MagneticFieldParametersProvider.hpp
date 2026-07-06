#pragma once

#include <memory>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TrackingPipeline/MagneticField/MagneticFieldStoreCollection.hpp"

namespace E320 {

class E320MagneticFieldParametersProvider {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Name of the input tree
    std::string treeName;
    /// The names of the input files
    std::vector<std::string> filePaths;
  };

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  explicit E320MagneticFieldParametersProvider(const Config& config);

  /// Get field store collection
  std::shared_ptr<MagneticFieldStoreCollection>
  getMagneticFieldStoreCollection();

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// The config class
  Config m_cfg;

  /// The input file and tree
  TFile* m_file = nullptr;
  TTree* m_tree = nullptr;
  TChain* m_chainOwner = nullptr;

 protected:
  /// Alignment store
  std::shared_ptr<MagneticFieldStoreCollection> m_storeCollection = nullptr;

  /// Quad gradients T/m
  double m_quad1Grad = 0;
  double m_quad2Grad = 0;
  double m_quad3Grad = 0;

  /// Quad poitions mm
  double m_quad1CenterPrimary = 0;
  double m_quad2CenterPrimary = 0;
  double m_quad3CenterPrimary = 0;

  double m_quad1CenterLong = 0;
  double m_quad2CenterLong = 0;
  double m_quad3CenterLong = 0;

  double m_quad1CenterShort = 0;
  double m_quad2CenterShort = 0;
  double m_quad3CenterShort = 0;

  /// XCOR strength T
  double m_xCorrectorStrength = 0;

  /// Dipole strength T
  double m_dipoleStrength = 0;

  /// Event number handle
  ULong64_t m_eventId = 0;
};

}  // namespace E320
