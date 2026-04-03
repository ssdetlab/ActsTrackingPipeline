#pragma once

#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <cstddef>
#include <utility>
#include <vector>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Simulation/IMeasurementGenerator.hpp"

class UniformBackgroundCreator : public IMeasurementGenerator {
 public:
  using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                      Acts::Experimental::DetectorNavigator>;

  using TrackParameters = Acts::CurvilinearTrackParameters;

  /// @brief Nested configuration struct
  struct Config {
    std::size_t nMeasurements;
    std::vector<const Acts::Surface*> surfaces;
    std::pair<double, double> resolution;
  };

  /// @brief Constructor
  explicit UniformBackgroundCreator(const Config& config);

  /// @brief Propagate track and create measurements
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
};
