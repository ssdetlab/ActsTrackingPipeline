#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/BinUtility.hpp"

#include "G4Transform3D.hh" 
#include "G4ThreeVector.hh"

#include <numeric>
#include <map>

/// Namespace for the E320 geometry
/// parameters
namespace E320Geometry {
    using namespace Acts::UnitLiterals;

    struct GeometryOptions {
        /// Constants for the translation 
        /// to the local coordinates
        const float chipSizeY = 29.94176_mm;
        const float chipSizeX = 13.76256_mm;

        // Chip Id to the translation in the x direction
        const std::map<std::int32_t, Acts::ActsScalar> 
        chipY{ 
            {0,  90.3_mm}, 
            {1, 120.4_mm}, 
            {2, 150.5_mm}, 
            {3, 180.6_mm}, 
            {4, 210.7_mm}, 
            {5, 240.8_mm}, 
            {6, 270.9_mm}, 
            {7, 301.0_mm}, 
            {8, 331.1_mm}
        };

        // Stave Id to the translation in the z direction
        const std::map<std::int32_t, Acts::ActsScalar>
        staveZ{
            {0, 16567.0_mm},
            {1, 16667.0_mm},
            {2, 16767.0_mm},
            {3, 16867.0_mm},
        };

        // All the staves are at the same y position
        Acts::ActsScalar chipX = -0.61872_mm;

        // Positions of the volumes
        // encapsulating the layers
        const std::vector<Acts::ActsScalar> layerZPositions{
            staveZ.at(0), 
            staveZ.at(1),
            staveZ.at(2),
            staveZ.at(3)
        };

        // Tracker volume encapsulating the 
        // whole detector
        const Acts::Vector3 trackerTranslation{
            0_mm, 0_mm, 8550_mm};

        const std::vector<Acts::ActsScalar> trackerBounds = {
            chipSizeX/2 + 1000_mm, 
            chipY.at(8) + chipSizeY/2 + 1000_mm,
            8600_mm
        };

        // Dipole volume encapsulating the
        // magnetic field
        const Acts::Vector3 dipoleTranslation{
            0_mm, 0_mm, 13140_mm};

        const std::vector<Acts::ActsScalar> dipoleBounds = {
            chipSizeX/2 + 1000_mm, 
            chipY.at(8) + chipSizeY/2 + 1000_mm,
            1000_mm
        };

        // Quadrupole volume encapsulating the
        // magnetic field
        const Acts::Vector3 quad1Translation{
            0_mm, 0_mm, 4182.49_mm};

        const std::vector<Acts::ActsScalar> quad1Bounds = {
            chipSizeX/2 + 1000_mm, 
            chipY.at(8) + chipSizeY/2 + 1000_mm,
            1000_mm
        };

        const Acts::Vector3 quad2Translation{
            0_mm, 0_mm, 6406.62_mm};

        const std::vector<Acts::ActsScalar> quad2Bounds = {
            chipSizeX/2 + 1000_mm, 
            chipY.at(8) + chipSizeY/2 + 1000_mm,
            1000_mm
        };

        const Acts::Vector3 quad3Translation{
            0_mm, 0_mm, 8631.05_mm};

        const std::vector<Acts::ActsScalar> quad3Bounds = {
            chipSizeX/2 + 1000_mm, 
            chipY.at(8) + chipSizeY/2 + 1000_mm,
            1000_mm
        };

        // Dipole and quadrupole field parameters
        const std::tuple<Acts::Vector2,Acts::Vector4,Acts::Vector4> 
        dipoleParams = {
            {-30.0_mm, 30.0_mm},
            {-165.0_mm, 165.0_mm, 7.7_mm, 7.7_mm},
            {-457.0_mm, 457.0_mm, 25.0_mm, 25.0_mm}
        };

        const Acts::Vector3 quadrupolesParams = {
            4.0_T / 1_m, 
            -7.0_T / 1_m,
            4.0_T / 1_m
        };

        // PDC window volume
        const Acts::Vector3 pdcWindowTranslation{
            0_mm, 0_mm, 16549.7_mm};

        const std::vector<Acts::ActsScalar> pdcWindowBounds = {
            chipSizeX/2 + 1000_mm, 
            chipY.at(8) + chipSizeY/2 + 1000_mm, 
            2_mm
        };

        // Arm volume encapsulating the layers
        const Acts::Vector3 armTranslation{
            0_mm, 0_mm, (staveZ.at(3) + staveZ.at(0))/2};

        /// Layer volumes encapsulating the staves
        const Acts::ActsScalar deltaZ = 
            (staveZ.at(1) - staveZ.at(0))/10 + 1_mm;

        const std::vector<Acts::ActsScalar> layerBounds = {
            chipSizeX/2 + 1000_mm, 
            chipY.at(8) + chipSizeY/2 + 1000_mm, 
            deltaZ
        };

        /// Global rotation of the world volume
        /// in the Acts format for volumes
        const Acts::RotationMatrix3 actsToWorldRotation = 
            Acts::AngleAxis3(M_PI_2,
                Acts::Vector3(1., 0., 0.)).toRotationMatrix();
        
        /// Global translation of the world volume
        /// in the Acts format for volumes
        const Acts::Translation3 actsToWorldTranslation = 
            Acts::Translation3(0, 0, 0);
        
        /// In case we'll need something less trivial 
        const Acts::Transform3 actsToWorld = 
            Acts::Transform3(actsToWorldRotation * actsToWorldTranslation);

        /// Global rotation of the world volume
        /// in the Geant4 format for surfaces
        const G4Transform3D g4ToWorld = G4Transform3D(
            CLHEP::HepRotationX(M_PI_2), 
            G4ThreeVector(0, 0, 0));

        // Material binning for the surfaces
        const Acts::BinUtility materialBinningX = 
            Acts::BinUtility(
                128, 
                -chipSizeY/2, 
                chipSizeY/2, 
                Acts::closed, 
                Acts::BinningValue::binX);

        const Acts::BinUtility materialBinningY = 
            Acts::BinUtility(
                64, 
                -chipSizeX/2, 
                chipSizeX/2, 
                Acts::closed, 
                Acts::BinningValue::binY);
    };

} // namespace LUXEGeometry
