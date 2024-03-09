#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"

#include <map>

/// Namespace for the LUXE geometry
/// parameters
namespace LUXEGeometry {
    using namespace Acts::UnitLiterals;
    
    /// Constants for the translation 
    /// to the local coordinates
    float chipSizeX = 29.94176_mm;
    float chipSizeY = 13.76256_mm;

    float chipTranslationY = 0.61872_mm - chipSizeY/2.0;
    std::map<int, float> chipTranslationXEven {
        {0, 67.73_mm  - chipSizeX/2.0},
        {1, 97.83_mm  - chipSizeX/2.0},
        {2, 127.93_mm - chipSizeX/2.0},
        {3, 158.03_mm - chipSizeX/2.0},
        {4, 188.13_mm - chipSizeX/2.0},
        {5, 218.23_mm - chipSizeX/2.0},
        {6, 248.33_mm - chipSizeX/2.0},
        {7, 278.43_mm - chipSizeX/2.0},
        {8, 308.53_mm - chipSizeX/2.0}
    };
    std::map<int, float> chipTranslationXOdd {
        {0, 298.53_mm - chipSizeX/2.0},
        {1, 328.63_mm - chipSizeX/2.0},
        {2, 358.73_mm - chipSizeX/2.0},
        {3, 388.83_mm - chipSizeX/2.0},
        {4, 418.93_mm - chipSizeX/2.0},
        {5, 449.03_mm - chipSizeX/2.0},
        {6, 479.13_mm - chipSizeX/2.0},
        {7, 509.23_mm - chipSizeX/2.0},
        {8, 539.33_mm - chipSizeX/2.0}
    };

    struct GeometryOptions {

        std::vector<Acts::ActsScalar> layerZPositions = {3962_mm,
                                                         4062_mm,
                                                         4162_mm,
                                                         4262_mm};
        // Can be set to zero if there are no IN and OUT staves for a single layer
        Acts::ActsScalar deltaZ = 12_mm;

        const std::vector<Acts::ActsScalar> detectorBounds = {100_mm,100_mm,5000_mm};

        const std::vector<Acts::ActsScalar> layerBounds = {100_mm, 100_mm, 14_mm};
    };

} // namespace LUXEGeometry