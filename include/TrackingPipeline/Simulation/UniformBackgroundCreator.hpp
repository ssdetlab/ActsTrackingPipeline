#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Simulation/IMeasurementGenerator.hpp"

/// @brief Class generating uniformly distributed background on sensitive surfaces
class UniformBackgroundCreator : public IMeasurementGenerator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Number of measurements to generate per surface
    std::size_t nMeasurementsPerSurface;
    /// Surfaces to generate on
    std::vector<const Acts::Surface*> surfaces;
    /// x-y resolution of the generated hits
    std::pair<double, double> resolution;
  };

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  explicit UniformBackgroundCreator(const Config& cfg);

  /// @brief Propagate track and create measurements
  ///
  /// @param ctx current algorithm context
  /// @param rng random engine for sampling
  /// @param id user-specified source id
  /// @param sourceLinks container to store created source links
  /// @param simClusters container to store created sim clusters
  /// @param sourceLinksIndices container to store created source link indices
  ///
  /// @return number of created measurements
  std::size_t gen(const AlgorithmContext& ctx, RandomEngine& rng,
                  std::size_t id, std::vector<Acts::SourceLink>& sourceLinks,
                  SimClusters& simClusters,
                  std::vector<std::size_t>& sourceLinksIndices) const override;

  /// @brief Readonly access to the config parameters
  const Config& config() const;

 private:
  /// Configuration
  Config m_cfg;

  /// Covariance
  Acts::SquareMatrix2 m_cov;

  /// Unfiorm distribution instance
  mutable std::uniform_real_distribution<double> m_uniform{0.0, 1.0};
};
