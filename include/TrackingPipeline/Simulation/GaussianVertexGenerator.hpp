#pragma once

#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/detail/NormalRandomVariable.hpp"

/// @brief Class generating gaussian distributed vertex vectors
/// with user-defined mean and covariance matrix
class GaussianVertexGenerator : public IVertexGenerator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Vertex mean
    Acts::Vector3 mean;
    /// Vertex covariance matrix
    Acts::SquareMatrix3 cov;
  };

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  explicit GaussianVertexGenerator(const Config& cfg);

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

  /// Normal distribution instance
  NormalRandomVariable m_normal;
};
