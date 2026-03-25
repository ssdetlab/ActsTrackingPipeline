#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include <Acts/Geometry/GeometryIdentifier.hpp>

#include <unordered_map>
#include <utility>

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

/// @brief Class that digitizes hits based on the provided
/// resolution assuming Gaussian smearing
class SurfaceRangedDigitizer : public IDigitizer {
 public:
  using Resolution = std::pair<double, double>;

  struct Config {
    std::unordered_map<Acts::GeometryIdentifier, Resolution> resolutions;
  };

  SurfaceRangedDigitizer(const Config& config);

  std::pair<Acts::Vector2, Acts::SquareMatrix2> genCluster(
      RandomEngine& rng, const Acts::GeometryIdentifier& geoId,
      const Acts::Vector2& pos) const override;

 private:
  Config m_cfg;

  std::unordered_map<Acts::GeometryIdentifier, Acts::SquareMatrix2> m_covs;
};
