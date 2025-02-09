#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

/// @brief Class that digitizes hits based on the provided
/// resolution assuming Gaussian smearing
struct SimpleDigitizer : public IDigitizer {
  std::pair<double, double> resolution;

  std::pair<Acts::SquareMatrix2, Acts::Vector2> genCluster(
      RandomEngine& rng, Acts::GeometryIdentifier /*geoId*/,
      Acts::Vector2 pos) const override {
    std::normal_distribution<double> normalDist(0., 1.);

    Acts::Vector2 stdDev = {resolution.first, resolution.second};
    Acts::Vector2 digLocal = pos + stdDev.cwiseProduct(Acts::Vector2(
                                       normalDist(rng), normalDist(rng)));

    Acts::SquareMatrix2 cov = stdDev.cwiseProduct(stdDev).asDiagonal();

    return {cov, digLocal};
  }
};
