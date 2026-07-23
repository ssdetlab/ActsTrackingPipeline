#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <unordered_map>
#include <utility>

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

/// @brief surface-specific Gaussian-smearing digitizer
class SurfaceRangedDigitizer : public IDigitizer {
 public:
  using Resolution = std::pair<double, double>;

  struct Config {
    std::unordered_map<std::size_t, Resolution> resolutions;
  };

  explicit SurfaceRangedDigitizer(const Config& config);

  std::pair<Acts::Vector2, Acts::SquareMatrix2> digitize(
      RandomEngine& rng, const Acts::GeometryIdentifier& geoId,
      const Acts::Vector2& pos) const override;

 private:
  Config m_cfg;

  std::unordered_map<std::size_t, Acts::SquareMatrix2> m_covs;

  mutable std::normal_distribution<double> m_normal{0., 1.};
};
