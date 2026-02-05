#pragma once

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Interface for generating backgorund hits
struct IMeasurementGenerator {
  virtual std::tuple<std::vector<Acts::SourceLink>, SimClusters> gen(
      const AlgorithmContext& ctx, RandomEngine& rng, std::size_t id) const = 0;
};
