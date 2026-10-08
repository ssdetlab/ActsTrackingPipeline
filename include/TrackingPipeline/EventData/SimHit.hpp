#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameters.hpp"

#include <vector>

/// @brief Track hit with truth information
struct SimHit {
  /// True parameters at the surface
  Acts::BoundVector truthParameters;
  /// Global hit position
  Acts::Vector3 globalPosition;
  /// True IP parameters
  Acts::CurvilinearTrackParameters ipParameters;
  /// True track Ids
  int trackId;
  /// True parent track Ids
  int parentTrackId;
  /// Run ID for unique identification
  int runId;
};

/// @brief Collection of SimHits
using SimHits = std::vector<SimHit>;
