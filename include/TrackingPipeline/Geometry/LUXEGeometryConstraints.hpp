#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/BinUtility.hpp"

#include <map>

#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

/// Namespace for the LUXE geometry
/// parameters
namespace LUXEGeometry {
using namespace Acts::UnitLiterals;

struct GeometryOptions {
  /// Constants for the translation
  /// to the local coordinates
  const float chipSizeX = 29.94176_mm;
  const float chipSizeY = 13.76256_mm;

  // Chip Id to the translation in the x direction
  // for the inner stave
  const std::map<std::int32_t, double> chipXEven{
      {0, 67.73_mm},  {1, 97.83_mm},  {2, 127.93_mm},
      {3, 158.03_mm}, {4, 188.13_mm}, {5, 218.23_mm},
      {6, 248.33_mm}, {7, 278.43_mm}, {8, 308.53_mm}};

  // Chip Id to the translation in the x direction
  // for the outer stave
  const std::map<std::int32_t, double> chipXOdd{
      {0, 298.53_mm}, {1, 328.63_mm}, {2, 358.73_mm},
      {3, 388.83_mm}, {4, 418.93_mm}, {5, 449.03_mm},
      {6, 479.13_mm}, {7, 509.23_mm}, {8, 539.33_mm}};

  // Stave Id to the translation in the z direction
  const std::map<std::int32_t, double> staveZ{
      {0, 3962.0125_mm}, {1, 3950.0125_mm}, {2, 4062.0125_mm},
      {3, 4050.0125_mm}, {4, 4162.0125_mm}, {5, 4150.0125_mm},
      {6, 4262.0125_mm}, {7, 4250.0125_mm}};

  // All the staves are at the same y position
  double chipY = 0.61872_mm;

  // Positions of the volumes
  // encapsulating the layers
  const std::vector<double> layerZPositions{
      (staveZ.at(0) + staveZ.at(1)) / 2, (staveZ.at(2) + staveZ.at(3)) / 2,
      (staveZ.at(4) + staveZ.at(5)) / 2, (staveZ.at(6) + staveZ.at(7)) / 2};

  // Tracker volume encapsulating the
  // whole detector
  const Acts::Vector3 trackerTranslation{0_mm, 0_mm, 2150_mm};

  const std::vector<double> trackerBounds = {
      chipXOdd.at(8) + chipSizeX / 2 + 100_mm, chipSizeY / 2 + 100_mm, 2200_mm};

  // Dipole volume encapsulating the
  // magnetic field
  const Acts::Vector3 dipoleTranslation{0_mm, 0_mm, 2050_mm};

  const std::vector<double> dipoleBounds = {
      chipXOdd.at(8) + chipSizeX / 2 + 100_mm, chipSizeY / 2 + 100_mm, 600_mm};

  const std::vector<double> constantFieldDelta = {0_mm, 0_mm, 120_mm};

  // Arm volume encapsulating the layers
  const Acts::Vector3 armTranslation{0_mm, 0_mm,
                                     (staveZ.at(7) + staveZ.at(0)) / 2};

  /// Layer volumes encapsulating the staves
  const double deltaZ = (staveZ.at(0) - staveZ.at(1)) / 2 + 1_mm;

  const std::vector<double> layerBounds = {
      chipXOdd.at(8) + chipSizeX / 2 + 100_mm, chipSizeY / 2 + 100_mm, deltaZ};

  /// Global rotation of the world volume
  /// in the Acts format for volumes
  const Acts::RotationMatrix3 actsToWorldRotation =
      Acts::AngleAxis3(M_PI_2, Acts::Vector3(1., 0., 0.)).toRotationMatrix();

  /// Global translation of the world volume
  /// in the Acts format for volumes
  const Acts::Translation3 actsToWorldTranslation = Acts::Translation3(0, 0, 0);

  /// In case we'll need something less trivial
  const Acts::Transform3 actsToWorld =
      Acts::Transform3(actsToWorldRotation * actsToWorldTranslation);

  /// Global rotation of the world volume
  /// in the Geant4 format for surfaces
  const G4Transform3D g4ToWorld =
      G4Transform3D(CLHEP::HepRotationX(M_PI_2), G4ThreeVector(0, 0, 0));

  // Path to the material file
  const std::string materialPath = "";

  // Material binning for the surfaces
  const Acts::BinUtility materialBinningX =
      Acts::BinUtility(16, -chipSizeX / 2, chipSizeX / 2, Acts::closed,
                       Acts::BinningValue::binX);

  const Acts::BinUtility materialBinningY = Acts::BinUtility(
      8, -chipSizeY / 2, chipSizeY / 2, Acts::closed, Acts::BinningValue::binY);
};

}  // namespace LUXEGeometry
