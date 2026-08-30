#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Interface for generating vertex vectors
struct IVertexGenerator {
  /// @brief Generate vertex vector
  ///
  /// @param rng random engine for sampling
  ///
  /// @return Generated momentum vector
  virtual Acts::Vector3 genVertex(RandomEngine& rng) const = 0;

  /// @brief Get vertex vector covariance
  ///
  /// @return Vertex covariance matrix in the format [x, y, z]
  virtual Acts::SquareMatrix3 getVertexCovariance() const = 0;

  /// @brief Get vertex vector mean
  ///
  /// @return Vertex mean
  virtual Acts::Vector3 getVertexMean() const = 0;
};
