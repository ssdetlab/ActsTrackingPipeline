#pragma once

#include <cstddef>
#include <unordered_map>

#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"

namespace E320 {

/// @brief E320-specific magnetic field store
struct E320MagneticFieldStore : public MagneticFieldStore {
  /// Quadrupole triplet focuing energy
  double quadTripletFocusingEnergy = 0;
  /// Quadrupole triplet object z
  double quadTripletObjectZ = 0;
  /// Quadrupole triplet image z
  double quadTripletImageZ = 0;
  /// Quadrupole triplet M12 transport matrix element
  double quadTripletM12 = 0;
  /// Quadrupole triplet M34 transport matrix element
  double quadTripletM34 = 0;
};

using E320MagneticFieldStores =
    std::unordered_map<std::size_t, std::shared_ptr<E320MagneticFieldStore>>;

}  // namespace E320
