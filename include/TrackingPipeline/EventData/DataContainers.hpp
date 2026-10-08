#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"

#include <cstddef>
#include <memory>
#include <vector>

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
