#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <stdexcept>

namespace detail {

inline Acts::BinningValue directionToBinningValue(const Acts::Vector3& dir) {
  if (dir == Acts::Vector3::UnitX()) {
    return Acts::BinningValue::binX;
  } else if (dir == Acts::Vector3::UnitY()) {
    return Acts::BinningValue::binY;
  } else if (dir == Acts::Vector3::UnitZ()) {
    return Acts::BinningValue::binZ;
  } else {
    throw std::runtime_error("Invalid vector to binning value conversion");
  }
}

inline Acts::Vector3 binningValueToDirection(
    const Acts::BinningValue& binningValue) {
  switch (binningValue) {
    case Acts::BinningValue::binX:
      return Acts::Vector3::UnitX();
    case Acts::BinningValue::binY:
      return Acts::Vector3::UnitY();
    case Acts::BinningValue::binZ:
      return Acts::Vector3::UnitZ();
    default:
      throw std::runtime_error("Invalid binning value to vector conversion");
  }
}

inline std::size_t binningValueToIndex(const Acts::BinningValue& binningValue) {
  switch (binningValue) {
    case Acts::BinningValue::binX:
      return 0;
    case Acts::BinningValue::binY:
      return 1;
    case Acts::BinningValue::binZ:
      return 2;
    default:
      throw std::runtime_error("Invalid binning value to index conversion");
  }
}

};  // namespace detail
