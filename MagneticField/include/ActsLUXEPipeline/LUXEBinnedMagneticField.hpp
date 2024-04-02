#pragma once


#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"

#include <vector>
#include <functional>
#include <tuple>

namespace LUXEMagneticField {

using vBins = std::vector<std::double_t>;
using eBins = std::tuple<std::double_t, std::double_t, std::size_t>;

using vAxis = Acts::detail::VariableAxis;
using eAxis = Acts::detail::EquidistantAxis;

using vGrid = Acts::Grid<Acts::Vector3, vAxis, vAxis, vAxis>;
using eGrid = Acts::Grid<Acts::Vector3, eAxis, eAxis, eAxis>;

/// Grid options for the 
/// variable axis
struct vGridOptions {
    vBins xBins;
    vBins yBins;
    vBins zBins;
};

/// Grid options for the 
/// equidistant axis
struct eGridOptions {
    eBins xBins;
    eBins yBins;
    eBins zBins;
};

/// Field and position transformations
using posTransform = std::function<
    Acts::Vector3(const Acts::Vector3)>;

using fieldTransform = std::function<
    Acts::Vector3(const Acts::Vector3,
        const Acts::Vector3)>;

using fieldVal = std::function<
    Acts::Vector3(const std::array<std::double_t, 3>)>;

/// @brief Construct binned magnetic field
/// with the variable axis
///
/// @param mFieldVal The magnetic field value
/// provider
/// @param transformPos The position transformation
/// @param transformBField The magnetic field 
/// transformation
/// @param gridOpt The grid options for the
/// variable axis
/// @param mctx The magnetic field context
/// @return an interpolated B field based
/// on a linear interpolation from values
/// calculated on a grid
Acts::InterpolatedBFieldMap<vGrid> buildBinnedBField(
    const Acts::MagneticFieldProvider& mFieldVal,
    const posTransform& transformPos,
    const fieldTransform& transformBField,
    const vGridOptions& gridOpt,
    const Acts::MagneticFieldContext& mctx);

/// @brief Construct binned magnetic field
/// with the equidistant axis
///
/// @param mFieldVal The magnetic field value
/// provider
/// @param transformPos The position transformation
/// @param transformBField The magnetic field 
/// transformation
/// @param gridOpt The grid options for the
/// equidistant axis
/// @param mctx The magnetic field context
/// @return an interpolated B field based
/// on a linear interpolation from values
/// calculated on a grid
Acts::InterpolatedBFieldMap<eGrid> buildBinnedBField(
    const Acts::MagneticFieldProvider& mFieldVal,
    const posTransform& transformPos,
    const fieldTransform& transformBField,
    const eGridOptions& gridOpt,
    const Acts::MagneticFieldContext& mctx); 

} // namespace LUXEMagneticField
