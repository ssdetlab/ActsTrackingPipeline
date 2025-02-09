#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackProxy.hpp"

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

/// @brief Track candidate to be passed to the KF
struct TrackCandidate {
  /// Source links related
  /// to the seed measurements
  std::vector<Acts::SourceLink> sourceLinks;
  /// IP parameters
  Acts::CurvilinearTrackParameters ipParameters;
  /// Track Id
  std::int32_t trackId;
  /// Predicted chi2 at the tracking layers
  std::vector<double> predictedChi2;
  /// Filtered chi2 at the tracking layers
  std::vector<double> filteredChi2;
};

/// @brief Collection of TrackCandidates
using TrackCandidates = std::vector<TrackCandidate>;

/// @brief A track with observables and uncertainties
template <typename container_t, typename trajectory_t>
struct Tracks {
  using TrackContainer =
      Acts::TrackContainer<container_t, trajectory_t, std::shared_ptr>;

  using IdContainer = std::vector<std::int32_t>;

  using ConstTrackProxy =
      Acts::TrackProxy<container_t, trajectory_t, std::shared_ptr, true>;

  TrackContainer tracks;
  IdContainer trackIds;

  std::pair<std::int32_t, ConstTrackProxy> getByIndex(std::int32_t i) {
    return {trackIds.at(i), tracks.getTrack(i)};
  }

  std::pair<std::int32_t, ConstTrackProxy> getByTrackId(std::int32_t i) {
    auto it = std::find(trackIds.begin(), trackIds.end(), i);
    if (it != trackIds.end()) {
      return {i, tracks.getTrack(std::distance(trackIds.begin(), it))};
    } else {
      return {i, nullptr};
    }
  }

  std::pair<std::int32_t, ConstTrackProxy> begin() {
    return {trackIds.at(0), tracks.getTrack(0)};
  }

  std::pair<std::int32_t, ConstTrackProxy> end() {
    return {trackIds.at(trackIds.size() - 1),
            tracks.getTrack(trackIds.size() - 1)};
  }

  std::size_t size() { return trackIds.size(); }
};
