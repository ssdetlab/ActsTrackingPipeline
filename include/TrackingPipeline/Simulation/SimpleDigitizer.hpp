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

  SimpleDigitizer(const Config& config);

  Acts::Vector2 genCluster(RandomEngine& rng,
                           const Acts::GeometryIdentifier& geoId,
                           const Acts::Vector2& pos) const override;

  Acts::SquareMatrix2 getCovariance(
      const Acts::GeometryIdentifier& geoId) const override;

 private:
  Config m_cfg;

  Acts::SquareMatrix2 m_cov;
};
