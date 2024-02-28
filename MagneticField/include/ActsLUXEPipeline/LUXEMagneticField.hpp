#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"


#include <iostream>

namespace LUXEMagneticField {

/// @brief Construct the test magnetic field
/// of a simple dipole located at the origin
///
/// @return an interpolated B field based
/// on a linear interpolation from values
/// calculated on a grid
using Grid_t =
        Acts::Grid<Acts::Vector3, Acts::detail::EquidistantAxis
                , Acts::detail::EquidistantAxis
                , Acts::detail::EquidistantAxis>;
using BField_t = Acts::InterpolatedBFieldMap<Grid_t>;

BField_t buildLUXEBField();

} // namespace LUXEGeometry
