#pragma once

#include <cstddef>

#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"

class MagneticFieldParametersContext {
 public:
  explicit MagneticFieldParametersContext(
      const MagneticFieldStoreCollection& collection)
      : m_collection(collection) {}

  MagneticFieldStore& magneticFieldStore(std::size_t eventNumber) {
    auto ub = m_collection.upper_bound(eventNumber);
    ub--;
    return *ub->second;
  }

  const MagneticFieldStore& magneticFieldStore(std::size_t eventNumber) const {
    auto ub = m_collection.upper_bound(eventNumber);
    ub--;
    return *ub->second;
  }

  std::shared_ptr<MagneticFieldStore> magneticFieldStorePtr(
      std::size_t eventNumber) {
    auto ub = m_collection.upper_bound(eventNumber);
    ub--;
    return ub->second;
  }

  const std::shared_ptr<MagneticFieldStore> magneticFieldStorePtr(
      std::size_t eventNumber) const {
    auto ub = m_collection.upper_bound(eventNumber);
    ub--;
    return ub->second;
  }

 private:
  MagneticFieldStoreCollection m_collection;
};
