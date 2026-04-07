#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"

#include <cstddef>
#include <memory>
#include <vector>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

///-----------------------------------------------
/// Sim data containers

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

/// @brief index-based seed container
struct IndexSeed {
  /// Source links indices related
  /// to the seed measurements
  std::vector<std::size_t> sourceLinkIndices;
  /// IP parameters index
  std::size_t originParametersIndex;
  /// Track Id
  int trackId;
};

/// @brief collection of Seeds
using IndexSeeds = std::vector<IndexSeed>;

/// @brief index-based track container
struct IndexTrack {
  std::size_t trackIndex;
  std::size_t originParametersGuessIndex;
  int trackId;
};

/// @brief collection of Tracks
using IndexTracks = std::vector<IndexTrack>;

/// -----------------------------------------------
/// Legacy data containers (slowly phased out)

/// @brief SourceLink-based seed container
struct Seed {
  /// Source links related
  /// to the seed measurements
  std::vector<Acts::SourceLink> sourceLinks;
  /// IP parameters
  Acts::CurvilinearTrackParameters ipParameters;
  /// Track Id
  int trackId;
};

/// @brief Collection of Seeds
using Seeds = std::vector<Seed>;

/// @brief Collection of Tracks
using ActsTracks =
    Acts::TrackContainer<Acts::VectorTrackContainer,
                         Acts::VectorMultiTrajectory, std::shared_ptr>;
struct Tracks {
  ActsTracks tracks;
  std::vector<int> trackIds;
  std::vector<Acts::CurvilinearTrackParameters> ipParametersGuesses;

  std::size_t size() const { return tracks.size(); }
};
