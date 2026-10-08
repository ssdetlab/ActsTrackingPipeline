#pragma once

#include "Acts/EventData/SourceLink.hpp"

#include <vector>

#include "TrackingPipeline/EventData/SimHit.hpp"

/// @brief Cluster with truth information
struct SimCluster {
  /// Observable parameters
  Acts::SourceLink sourceLink;
  /// Truth parameters
  SimHits truthHits;
  /// Is Signal flag
  bool isSignal;
};

/// @brief Collection of SimClusters
using SimClusters = std::vector<SimCluster>;
