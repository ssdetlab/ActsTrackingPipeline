#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"

#include "CLHEP/Vector/Rotation.h" 

#include <numeric>
#include <map>

/// Namespace for the LUXE geometry
/// parameters
namespace LUXEGeometry {
    using namespace Acts::UnitLiterals;

    struct GeometryOptions {
        /// Constants for the translation 
        /// to the local coordinates
        const float chipSizeX = 29.94176_mm;
        const float chipSizeY = 13.76256_mm;

        const float chipTranslationY = 0.61872_mm - chipSizeY/2.0;
        const std::map<int, float> chipTranslationXEven {
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
        const std::map<int, float> chipTranslationXOdd {
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

        const std::vector<Acts::ActsScalar> layerZPositions = 
            {3956.0125_mm, 4056.0125_mm, 4156.0125_mm, 4256.0125_mm};

        const Acts::Vector3 postironArmTranslation{
            305_mm, 0_mm, std::accumulate(layerZPositions.begin(), 
                layerZPositions.end(), 0_mm)/layerZPositions.size()};

        /// Can be set to zero if there are 
        /// no IN and OUT staves for a single layer
        const Acts::ActsScalar deltaZ = 7_mm;

        const std::vector<Acts::ActsScalar> postironArmBounds = 
            {255_mm, 10_mm, 160_mm};

        const std::vector<Acts::ActsScalar> layerBounds = 
            {255_mm, 10_mm, deltaZ};

        const std::vector<std::pair<Acts::ActsScalar,Acts::ActsScalar>> MagneticFieldBounds =
                {std::make_pair(-1000_mm,1000_mm),
                 std::make_pair(1450_mm,2650_mm),
                 std::make_pair(-100_mm,100_mm)};

        /// Global rotation of the world volume
        /// in the Acts format for volumes
        const Acts::RotationMatrix3 actsWorldRotation = 
            Acts::AngleAxis3(-M_PI_2,
                Acts::Vector3(1., 0., 0.)).toRotationMatrix();
        
        /// Global rotation of the world volume
        /// in the Geant4 format for surfaces
        const CLHEP::HepRotationX g4WorldRotation =
            CLHEP::HepRotationX(-M_PI_2);
    };

} // namespace LUXEGeometry
