#pragma once

#include <cstddef>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TRandom.h"
#include "TrackingPipeline/Simulation/IMultiplicityGenerator.hpp"

/// @brief Class generating multiplicity sampled from a ROOT TH1I hist
class RootHistMultiplicityGenerator : public IMultiplicityGenerator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Path to input file
    std::string filePath;
    /// Hist name
    std::string histName;
  };

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  explicit RootHistMultiplicityGenerator(const Config& cfg);

  /// @brief Destructor
  ~RootHistMultiplicityGenerator();

  /// @brief Generate event multiplicity
  ///
  /// @param rng random engine for sampling
  ///
  /// @return Number of tracks/measurement to generate
  std::size_t genMultiplicity(RandomEngine& rng) const override;

  /// @brief Get multiplicity standard error
  ///
  /// @return Multiplicity standard error
  double getMultiplicityStdDev() const override;

  /// @brief Get multiplicity mean
  ///
  /// @return Multiplicity mean
  double getMultiplicityMean() const override;

 private:
  /// Configuration
  Config m_cfg;

  /// File handle
  TFile* m_file = nullptr;

  /// Sampling hist
  TH1I* m_hist = nullptr;

  /// ROOT rng instance
  TRandom* m_rng = nullptr;
};
