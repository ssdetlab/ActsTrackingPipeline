#pragma once

#include <cstddef>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Interface for generating backgorund hits
struct IMeasurementGenerator {
  virtual std::size_t gen(
      const AlgorithmContext& ctx, RandomEngine& rng, std::size_t id,
      std::vector<Acts::SourceLink>& sourceLinks, SimClusters& simClusters,
      std::vector<std::size_t>& sourceLinksIndices) const = 0;
};
