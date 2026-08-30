#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

/// @brief Class generating uniformly distributed vertex vectors
/// in user-defined 3D bounds
class UniformVertexGenerator : public IVertexGenerator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Lower vertex bounds
    Acts::Vector3 mins;
    /// Higher vertex bounds
    Acts::Vector3 maxs;
  };

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  explicit UniformVertexGenerator(const Config& cfg);

  /// @brief Generate vertex vector
  ///
  /// @param rng random engine for sampling
  ///
  /// @return Generated momentum vector
  Acts::Vector3 genVertex(RandomEngine& rng) const override;

  /// @brief Get vertex vector covariance
  ///
  /// @return Vertex covariance matrix in the format [x, y, z]
  Acts::SquareMatrix3 getVertexCovariance() const override;

  /// @brief Get vertex vector mean
  ///
  /// @return Vertex mean
  Acts::Vector3 getVertexMean() const override;

 private:
  /// Configuration
  Config m_cfg;

  /// Vertex covariance
  Acts::SquareMatrix3 m_cov;

  /// Vertex mean
  Acts::Vector3 m_mean;

  /// Uniform distribution instance
  mutable std::uniform_real_distribution<double> m_uniform{0, 1};
};
