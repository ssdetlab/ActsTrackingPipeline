#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <cstddef>
#include <unordered_map>

/// @brief Struct wrapping the magnetic field store
/// and providing the inheritance base
struct MagneticFieldStore {
  /// Virtual destructor for inheritance interfacing
  virtual ~MagneticFieldStore() = default;

  /// Store map [field id, field cache]
  std::unordered_map<std::size_t, Acts::MagneticFieldProvider::Cache> store;
};

/// @brief Shorthand for a collection of stores
using MagneticFieldStores =
    std::unordered_map<std::size_t, std::shared_ptr<MagneticFieldStore>>;
