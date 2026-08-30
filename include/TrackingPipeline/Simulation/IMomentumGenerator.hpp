#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Interface for generating momentum vectors
struct IMomentumGenerator {
  /// @brief Generate momentum vector
  ///
  /// @param rng random engine for sampling
  ///
  /// @return Generated momentum vector
  virtual Acts::Vector3 genMomentum(RandomEngine& rng) const = 0;

  /// @brief Get momentum vector covariance
  ///
  /// @return Momentum covariance matrix in the format [dirX, dirY, dirZ, P]
  virtual Acts::SquareMatrix4 getMomentumCovariance() const = 0;

  /// @brief Get momentum vector mean
  ///
  /// @return Momentum mean
  virtual Acts::Vector3 getMomentumMean() const = 0;
};
