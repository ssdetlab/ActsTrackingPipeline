#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <unordered_map>

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

/// @brief Class that digitizes hits based on the simulated
/// cluster size distribution
class ClusterSizeBasedDigitizer : public IDigitizer {
 public:
  using SizeParameters = std::tuple<double, double, double>;

  struct Config {
    std::unordered_map<std::size_t, SizeParameters> clSizeProbsStdDevs;
  };

  ClusterSizeBasedDigitizer(const Config& config);

  std::pair<Acts::Vector2, Acts::SquareMatrix2> genCluster(
      RandomEngine& rng, const Acts::GeometryIdentifier& geoId,
      const Acts::Vector2& pos) const override;

 private:
  Config m_cfg;

  std::vector<double> m_clSizeProbs;

  Acts::SquareMatrix2 m_cov;
};
