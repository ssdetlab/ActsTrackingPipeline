#pragma once

#include <cstddef>
#include <vector>

/// @brief index-based track container
struct IndexTrack {
  /// Index inside the acts track container
  std::size_t trackIndex;
  /// Index inside the guess origin parameters container
  std::size_t originParametersGuessIndex;
  /// Track ID
  int trackId;
};

/// @brief collection of Tracks
using IndexTracks = std::vector<IndexTrack>;
