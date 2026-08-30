#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"

/// @brief Class generating momentum vectors with uniform density
/// magnitude and user-specified direction
class RangedUniformMomentumGenerator : public IMomentumGenerator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Momentum magnitude range
    std::vector<std::pair<double, double>> pRanges;
    /// Momentum direction
    Acts::Vector3 direction;
  };

  /// @brief Constructor
  ///
  /// @param config configuration struct
  explicit RangedUniformMomentumGenerator(const Config& config);

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

  /// Momentum mean
  Acts::Vector3 m_mean;

  /// Uniform int distribution instance
  mutable std::uniform_int_distribution<int> m_rangeSelect;
  
  /// Uniform real distribution instance
  mutable std::uniform_real_distribution<double> m_uniform{0, 1};
};
