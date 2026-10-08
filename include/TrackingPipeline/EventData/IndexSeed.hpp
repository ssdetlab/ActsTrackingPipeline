#pragma once

#include <cstddef>
#include <vector>

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
