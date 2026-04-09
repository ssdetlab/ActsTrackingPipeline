#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <cstddef>
#include <map>

struct MagneticFieldStore {
  std::unordered_map<std::size_t, Acts::MagneticFieldProvider::Cache> store;
};

using MagneticFieldStoreCollection =
    std::map<std::size_t, std::shared_ptr<MagneticFieldStore>>;
