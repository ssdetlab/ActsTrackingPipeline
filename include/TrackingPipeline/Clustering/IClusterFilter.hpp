#pragma once

#include "TrackingPipeline/EventData/DataContainers.hpp"

class IClusterFilter {
 public:
  virtual ~IClusterFilter() = default;

  virtual bool operator()(const Acts::GeometryContext& gctx,
                          const SimCluster& cluster) const = 0;
};
