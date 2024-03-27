#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "ActsLUXEPipeline/IMagneticField.hpp"
#include "ActsLUXEPipeline/MagneticFields.hpp"

#include <iostream>
#include <vector>
#include <functional>
#include <variant>
#include <tuple>
#include <type_traits>

namespace LUXEMagneticField {
using namespace Acts::UnitLiterals;
using bin_t = std::vector<double>;
using axis_t = Acts::detail::VariableAxis;

struct GridOptions {
    bin_t xBins;
    bin_t yBins;
    bin_t zBins;
};

/// 3D Variable binning
using Grid_t = Acts::Grid<Acts::Vector3, axis_t, axis_t, axis_t>;

using BField_t = Acts::InterpolatedBFieldMap<Grid_t>;

/// If any coordinate transformations are required
using transformationPos_t = std::function
        <Acts::Vector3(const Acts::Vector3)>;

using transformationBField_t = std::function
        <Acts::Vector3(const Acts::Vector3,
                       const Acts::Vector3)>;

using bFieldValue_t = std::function
        <Acts::Vector3(const std::array<double, 3>)>;

/// @brief Construct the test magnetic field
/// of a simple dipole located at the origin
///
/// @return an interpolated B field based
/// on a linear interpolation from values
/// calculated on a grid
BField_t buildLUXEBField(const transformationPos_t& transformPos,
                         const transformationBField_t& transformBField,
                         const GridOptions gridOpt,
                         const IMagneticField& magneticField);

} // namespace LUXEGeometry
