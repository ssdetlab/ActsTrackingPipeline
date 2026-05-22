#pragma once

#include <cstddef>

#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"

class MagneticFieldStoreCollection {
 public:
  explicit MagneticFieldStoreCollection(const MagneticFieldStores& collection)
      : m_stores(collection) {}

  MagneticFieldStore& magneticFieldStore(std::size_t id) {
    return *m_stores.at(id);
  }

  const MagneticFieldStore& magneticFieldStore(std::size_t id) const {
    return *m_stores.at(id);
  }

  std::shared_ptr<MagneticFieldStore> magneticFieldStorePtr(std::size_t id) {
    return m_stores.at(id);
  }

  const std::shared_ptr<MagneticFieldStore> magneticFieldStorePtr(
      std::size_t id) const {
    return m_stores.at(id);
  }

 private:
  MagneticFieldStores m_stores;
};
