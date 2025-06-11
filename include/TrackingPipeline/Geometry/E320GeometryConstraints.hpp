#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/BinUtility.hpp"

#include <cmath>
#include <map>

#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

/// Namespace for the E320 geometry
/// parameters
namespace E320Geometry {
using namespace Acts::UnitLiterals;

struct GeometryOptions {
  /// Constants for the translation
  /// to the local coordinates
  const float chipSizeY = 29.94176_mm;
  const float chipSizeX = 13.76256_mm;

  const float pixelSizeX = 26.88_um;
  const float pixelSizeY = 29.24_um;

  // Chip Id to the translation in the x direction
  double chipY = 93.775_mm;

  // Stave Id to the translation in the z direction
  const std::map<int, double> staveZ{{8, 16674.4_mm},
                                     {6, 16694.4_mm},
                                     {4, 16714.4_mm},
                                     {2, 16734.4_mm},
                                     {0, 16754.4_mm}};

  // All the staves are at the same y position
  double chipX = 0.61872_mm;

  // Positions of the volumes
  // encapsulating the layers
  const std::vector<double> layerZPositions{
      staveZ.at(8), staveZ.at(6), staveZ.at(4), staveZ.at(2), staveZ.at(0)};

  // Tracker volume encapsulating the
  // whole detector
  const Acts::Vector3 trackerTranslation{0_mm, 0_mm, 8550_mm};

  const std::vector<double> trackerBounds = {
      chipSizeX / 2 + 2000_mm, chipY + chipSizeY / 2 + 2000_mm, 8600_mm};

  // Be window poisiton
  const Acts::Vector3 beWindowTranslation{0_mm, 0_mm, -842_mm};
  const std::vector<double> beWindowBounds = {
      chipSizeX / 2 + 2000_mm, chipY + chipSizeY / 2 + 2000_mm, 10_mm};

  // Quadrupole volume encapsulating the
  // magnetic field
  const Acts::Vector3 quad1Translation{0_mm, 0_mm, 4160_mm};

  const std::vector<double> quad1Bounds = {
      chipSizeX / 2 + 2000_mm, chipY + chipSizeY / 2 + 2000_mm, 486.664_mm};

  const Acts::Vector3 quad2Translation{0_mm, 0_mm, 6390_mm};

  const std::vector<double> quad2Bounds = {
      chipSizeX / 2 + 2000_mm, chipY + chipSizeY / 2 + 2000_mm, 486.664_mm};

  const Acts::Vector3 quad3Translation{0_mm, 0_mm, 8610_mm};

  const std::vector<double> quad3Bounds = {
      chipSizeX / 2 + 2000_mm, chipY + chipSizeY / 2 + 2000_mm, 486.664_mm};

  // Dipole and quadrupole field parameters
  const std::tuple<Acts::Vector2, Acts::Vector4, Acts::Vector4> dipoleParams = {
      {-30.0_mm, 30.0_mm},
      {-165.0_mm, 165.0_mm, 7.7_mm, 7.7_mm},
      {-457.0_mm, 457.0_mm, 25.0_mm, 25.0_mm}};

  // Run 489
  /*const Acts::Vector3 quadrupolesParams = {-3.218_T / 1_m, 4.719_T / 1_m,*/
  /*                                         -3.218_T / 1_m};*/
  // Run 502
  const Acts::Vector3 quadrupolesParams = {-0.7637_T / 1_m, 2.855_T / 1_m,
                                           -0.7637_T / 1_m};

  // X-corrector volume encapsulating the
  // magnetic field
  const Acts::Vector3 xCorrectorTranslation{0_mm, 0_mm, 9994.63_mm};

  const std::vector<double> xCorrectorBounds = {
      chipSizeX / 2 + 2000_mm, chipY + chipSizeY / 2 + 2000_mm, 116.84_mm};

  // Dipole volume encapsulating the
  // magnetic field
  const Acts::Vector3 dipoleTranslation{0_mm, 0_mm, 13060.6_mm};

  const std::vector<double> dipoleBounds = {
      chipSizeX / 2 + 2000_mm, chipY + chipSizeY / 2 + 2000_mm, 1000_mm};

  // PDC window volume
  // TODO: Add machinery to move the PDC window
  const Acts::Vector3 pdcWindowTranslation{0_mm, 0_mm, 16549.7_mm};

  const std::vector<double> pdcWindowBounds = {
      chipSizeX / 2 + 2000_mm, chipY + chipSizeY / 2 + 2000_mm, 2_mm};

  // Arm volume encapsulating the layers
  const Acts::Vector3 armTranslation{0_mm, 0_mm,
                                     (staveZ.at(0) + staveZ.at(8)) / 2};

  /// Layer volumes encapsulating the staves
  const double deltaZ = (staveZ.at(6) - staveZ.at(8)) / 10 + 1_mm;

  const std::vector<double> layerBounds = {
      chipSizeX / 2 + 2000_mm, chipY + chipSizeY / 2 + 2000_mm, deltaZ};

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

  // Material binning for the surfaces
  const Acts::BinUtility materialBinningX =
      Acts::BinUtility(256, -chipSizeY / 2, chipSizeY / 2, Acts::closed,
                       Acts::BinningValue::binX);

  const Acts::BinUtility materialBinningY =
      Acts::BinUtility(128, -chipSizeX / 2, chipSizeX / 2, Acts::closed,
                       Acts::BinningValue::binY);
};

}  // namespace E320Geometry
