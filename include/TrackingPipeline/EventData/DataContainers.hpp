#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include <Acts/Utilities/Holders.hpp>

#include <memory>
#include <vector>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

///-----------------------------------------------
/// Simulated and truth data containers

struct SimHit {
  /// Source link for compatibility
  /// with some algorithms
  Acts::SourceLink sourceLink;
  /// Truth parameters
  Acts::BoundVector truthParameters;
  /// True IP parameters
  Acts::CurvilinearTrackParameters ipParameters;
  /// True track Ids
  std::int32_t trackId;
  /// True parent track Ids
  std::int32_t parentTrackId;
  /// Run ID for unique identification
  std::int32_t runId;
};

/// @brief Collection of SimHits
using SimHits = std::vector<SimHit>;

/// @brief Cluster with truth information
struct SimCluster {
  /// Observable parameters
  SimpleSourceLink sourceLink;
  /// Truth parameters
  SimHits truthHits;
  /// Is Signal flag
  bool isSignal;
};

/// @brief Collection of SimClusters
using SimClusters = std::vector<SimCluster>;

///-----------------------------------------------
/// Obserbable data containers

/// @brief Seed to be passed to the CKF
struct Seed {
  /// Source links related
  /// to the seed measurements
  std::vector<Acts::SourceLink> sourceLinks;
  /// IP parameters
  Acts::CurvilinearTrackParameters ipParameters;
  /// Track Id
  std::int32_t trackId;
};

/// @brief Collection of Seeds
using Seeds = std::vector<Seed>;

/// @brief Collection of Tracks
using Tracks =
    Acts::TrackContainer<Acts::VectorTrackContainer,
                         Acts::VectorMultiTrajectory, std::shared_ptr>;
