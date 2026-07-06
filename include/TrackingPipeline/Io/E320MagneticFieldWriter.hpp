#pragma once

#include <Acts/Utilities/Logger.hpp>

#include <memory>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TrackingPipeline/Infrastructure/IWriter.hpp"

namespace E320 {

class E320MagneticFieldWriter : public IWriter {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Name of the output tree
    std::string treeName;
    /// The names of the output file
    std::string filePath;
  };

  E320MagneticFieldWriter(const E320MagneticFieldWriter &) = delete;
  E320MagneticFieldWriter(const E320MagneticFieldWriter &&) = delete;

  /// @brief constructor
  E320MagneticFieldWriter(const Config &config, Acts::Logging::Level level);

  /// @brief finalize method
  ProcessCode finalize() override;

  /// @brief get the writer name
  std::string name() const override { return "E320MagneticFieldWriter"; }

  /// @brief write out data to the file
  ProcessCode write(const AlgorithmContext &ctx) override;

 protected:
  /// @brief access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// The config class
  Config m_cfg;

  /// The output file and tree
  TFile *m_file = nullptr;
  TTree *m_tree = nullptr;

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

 protected:
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

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};

}  // namespace E320
