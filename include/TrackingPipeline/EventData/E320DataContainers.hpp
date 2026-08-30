#pragma once

#include <cstddef>
#include <vector>

namespace E320 {

///-----------------------------------------------
/// Obserbable data containers

/// @brief index-based seed container
struct E320IndexSeed {
  /// Source links indices related
  /// to the seed measurements
  std::vector<std::size_t> sourceLinkIndices;
  /// IP parameters index
  std::size_t originParametersIndex;
  /// Track ID
  int trackId;
  /// Number of HT cell intersection counts
  std::size_t htXCount;
};

/// @brief collection of Seeds
using E320IndexSeeds = std::vector<E320IndexSeed>;

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
