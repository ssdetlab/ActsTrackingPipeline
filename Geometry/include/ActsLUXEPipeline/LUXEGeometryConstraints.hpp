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

        const std::map<std::int32_t, Acts::ActsScalar> 
        chipXEven{ 
            {0, 67.73_mm}, 
            {1, 97.83_mm}, 
            {2, 127.93_mm}, 
            {3, 158.03_mm}, 
            {4, 188.13_mm}, 
            {5, 218.23_mm}, 
            {6, 248.33_mm}, 
            {7, 278.43_mm}, 
            {8, 308.53_mm}
        };

        const std::map<std::int32_t, Acts::ActsScalar> 
        chipXOdd{
            {0, 298.53_mm}, 
            {1, 328.63_mm}, 
            {2, 358.73_mm}, 
            {3, 388.83_mm}, 
            {4, 418.93_mm}, 
            {5, 449.03_mm}, 
            {6, 479.13_mm}, 
            {7, 509.23_mm}, 
            {8, 539.33_mm}
        };

        const std::map<std::int32_t, Acts::ActsScalar>
        layerZ{
            {0, 3962.0125_mm},
            {1, 3950.0125_mm},
            {2, 4062.0125_mm},
            {3, 4050.0125_mm},
            {4, 4162.0125_mm},
            {5, 4150.0125_mm},
            {6, 4262.0125_mm},
            {7, 4250.0125_mm}
        };

        Acts::ActsScalar chipY = 0.61872_mm;

        const float chipTranslationY = chipY - chipSizeY/2.0;
        const std::map<int, float> chipTranslationXEven {
            {0, chipXEven.at(0) - chipSizeX/2.0},
            {1, chipXEven.at(1) - chipSizeX/2.0},
            {2, chipXEven.at(2) - chipSizeX/2.0},
            {3, chipXEven.at(3) - chipSizeX/2.0},
            {4, chipXEven.at(4) - chipSizeX/2.0},
            {5, chipXEven.at(5) - chipSizeX/2.0},
            {6, chipXEven.at(6) - chipSizeX/2.0},
            {7, chipXEven.at(7) - chipSizeX/2.0},
            {8, chipXEven.at(8) - chipSizeX/2.0}
        };

        const std::map<int, float> chipTranslationXOdd {
            {0, chipXOdd.at(0) - chipSizeX/2.0},
            {1, chipXOdd.at(1) - chipSizeX/2.0},
            {2, chipXOdd.at(2) - chipSizeX/2.0},
            {3, chipXOdd.at(3) - chipSizeX/2.0},
            {4, chipXOdd.at(4) - chipSizeX/2.0},
            {5, chipXOdd.at(5) - chipSizeX/2.0},
            {6, chipXOdd.at(6) - chipSizeX/2.0},
            {7, chipXOdd.at(7) - chipSizeX/2.0},
            {8, chipXOdd.at(8) - chipSizeX/2.0}
        };

        const std::vector<Acts::ActsScalar> layerZPositions{
            (layerZ.at(0) + layerZ.at(1))/2, 
            (layerZ.at(2) + layerZ.at(3))/2,
            (layerZ.at(4) + layerZ.at(5))/2,
            (layerZ.at(6) + layerZ.at(7))/2
        };

        const Acts::Vector3 postironArmTranslation{
            (chipXOdd.at(8) + chipXEven.at(0))/2, 
            0_mm, 
            2575_mm};

        const Acts::Vector3 electronArmTranslation{
            -(chipXOdd.at(8) + chipXEven.at(0))/2, 0_mm, 2575_mm};

        const Acts::Vector3 trackerTranslation{
            0_mm, 0_mm, 2575_mm};

        /// Can be set to zero if there are 
        /// no IN and OUT staves for a single layer
        const Acts::ActsScalar deltaZ = 
            (layerZ.at(0) - layerZ.at(1))/2 + 1_mm;

        const std::vector<Acts::ActsScalar> trackerBounds = 
            {(chipXOdd.at(8) - chipXEven.at(0))/2 + chipSizeX/2 + 1_mm, 
                chipSizeY/2 + 1_mm, 1725_mm};

        const std::vector<Acts::ActsScalar> armBounds = 
            {(chipXOdd.at(8) - chipXEven.at(0))/2 + chipSizeX/2 + 1_mm, 
                chipSizeY/2 + 1_mm, 1725_mm};

        const std::vector<Acts::ActsScalar> layerBounds = 
            {(chipXOdd.at(8) - chipXEven.at(0))/2 + 1_mm, 
                chipSizeY/2 + 1_mm, deltaZ};

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
