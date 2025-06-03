#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

/// @brief Class that digitizes hits based on the provided
/// resolution assuming Gaussian smearing
class SimpleDigitizer : public IDigitizer {
 public:
  struct Config {
    std::pair<double, double> resolution;
  };

  SimpleDigitizer(const Config& config) : m_cfg(config) {
    Acts::Vector2 stdDev = {m_cfg.resolution.first, m_cfg.resolution.second};
    m_cov = stdDev.cwiseProduct(stdDev).asDiagonal();
  }

  std::pair<Acts::SquareMatrix2, Acts::Vector2> genCluster(
      RandomEngine& rng, Acts::GeometryIdentifier /*geoId*/,
      Acts::Vector2 pos) const override {
    std::normal_distribution<double> normalDist(0., 1.);

    Acts::Vector2 stdDev = {m_cfg.resolution.first, m_cfg.resolution.second};
    Acts::Vector2 digLocal = pos + stdDev.cwiseProduct(Acts::Vector2(
                                       normalDist(rng), normalDist(rng)));
    return {m_cov, digLocal};
  }

 private:
  /// Configuration
  Config m_cfg;

  /// Covariance
  Acts::SquareMatrix2 m_cov;
};
