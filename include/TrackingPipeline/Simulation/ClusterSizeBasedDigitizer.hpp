#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <unordered_map>

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

/// @brief Class digitizing hits based on the simulated
/// cluster size distribution
class ClusterSizeBasedDigitizer : public IDigitizer {
 public:
  /// [Cluster-size sampling prob, resolution X, resolution Y] shorthand
  using SizeParameters = std::tuple<double, double, double>;

  /// @bried Nested config struct
  struct Config {
    /// Map of cluster sizes to their size parameters
    std::unordered_map<std::size_t, SizeParameters> clSizeProbsStdDevs;
  };

  /// @brief Constructor
  ///
  /// @param config configuration struct
  explicit ClusterSizeBasedDigitizer(const Config& config);

  /// @brief digitize measurement
  ///
  /// @param rng random engine for sampling
  /// @param geoId geometry ID of the measurement surface
  /// @param pos 2D measurement
  ///
  /// @return digitized measurement
  std::pair<Acts::Vector2, Acts::SquareMatrix2> digitize(
      RandomEngine& rng, const Acts::GeometryIdentifier& geoId,
      const Acts::Vector2& pos) const override;

 private:
  /// Configuration
  Config m_cfg;

  /// Cluster sizes probabilities of occurence
  std::vector<double> m_clSizeProbs;

  /// Measurement covariance
  Acts::SquareMatrix2 m_cov;

  /// Normal distribution instance
  mutable std::normal_distribution<double> m_normal{0., 1.};

  /// Dicrete distribution instance
  mutable std::discrete_distribution<std::size_t> m_discrete;
};
