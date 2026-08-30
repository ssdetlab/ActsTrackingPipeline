#pragma once

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"

/// @brief Class generating constant user-defined momentum vector
class StationaryMomentumGenerator : public IMomentumGenerator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Momentum vector to return
    Acts::Vector3 momentum;
  };

  /// @brief Constructor
  ///
  /// @param config configuration struct
  explicit StationaryMomentumGenerator(const Config& config);

  /// @brief Generate momentum vector
  ///
  /// @param rng random engine for sampling
  ///
  /// @return Generated momentum vector
  Acts::Vector3 genMomentum(RandomEngine& rng) const override;

  /// @brief Get momentum vector covariance
  ///
  /// @return Momentum covariance matrix in the format [dirX, dirY, dirZ, P]
  Acts::SquareMatrix4 getMomentumCovariance() const override;

  /// @brief Get momentum vector mean
  ///
  /// @return Momentum mean
  Acts::Vector3 getMomentumMean() const override;

 private:
  /// Configuration
  Config m_cfg;

  /// Momentum covariance matrix
  Acts::SquareMatrix4 m_cov;
};
