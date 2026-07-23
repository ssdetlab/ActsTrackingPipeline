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

  explicit SimpleDigitizer(const Config& config);

  std::pair<Acts::Vector2, Acts::SquareMatrix2> digitize(
      RandomEngine& rng, const Acts::GeometryIdentifier& geoId,
      const Acts::Vector2& pos) const override;

 private:
  Config m_cfg;

  Acts::SquareMatrix2 m_cov;

  mutable std::normal_distribution<double> m_normal{0., 1.};
};
