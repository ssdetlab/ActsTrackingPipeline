#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

/// @brief Class generating constant user-defined vertex vectors
class StationaryVertexGenerator : public IVertexGenerator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Vertex vector to generate
    Acts::Vector3 vertex;
  };

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  explicit StationaryVertexGenerator(const Config& cfg);

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
};
