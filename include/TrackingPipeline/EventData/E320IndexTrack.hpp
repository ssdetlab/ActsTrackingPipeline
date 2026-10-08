#pragma once

#include <cstddef>
#include <vector>

namespace E320 {

/// @brief index-based track container
struct E320IndexTrack {
  /// Index inside the acts track container
  std::size_t trackIndex;
  /// Index inside the guess origin parameters container
  std::size_t originParametersGuessIndex;
  /// Track ID
  int trackId;
  /// Number of HT cell intersection counts
  std::size_t xCount;
};

/// @brief collection of Tracks
using E320IndexTracks = std::vector<E320IndexTrack>;

}  // namespace E320
