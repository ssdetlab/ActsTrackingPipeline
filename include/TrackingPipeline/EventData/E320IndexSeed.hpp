#pragma once

#include <cstddef>
#include <vector>

namespace E320 {

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

}  // namespace E320
