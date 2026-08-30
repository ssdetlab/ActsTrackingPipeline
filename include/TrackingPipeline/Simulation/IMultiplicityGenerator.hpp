#pragma once

#include <cstddef>

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Interface for generating event multiplicity
struct IMultiplicityGenerator {
  /// @brief Generate event multiplicity
  ///
  /// @param rng random engine for sampling
  ///
  /// @return Number of tracks/measurement to generate
  virtual std::size_t genMultiplicity(RandomEngine& rng) const = 0;

  /// @brief Get multiplicity standard error
  ///
  /// @return Multiplicity standard error
  virtual double getMultiplicityStdDev() const = 0;

  /// @brief Get multiplicity mean
  ///
  /// @return Multiplicity mean
  virtual double getMultiplicityMean() const = 0;
};
