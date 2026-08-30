#pragma once

#include <cstddef>

#include "TrackingPipeline/Simulation/IMultiplicityGenerator.hpp"

/// @brief Class generating constant user-defined event multiplicity
class ConstantMultiplicityGenerator : public IMultiplicityGenerator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Event multiplicity
    std::size_t eventMultiplicity;
  };

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  explicit ConstantMultiplicityGenerator(const Config& cfg);

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
};
