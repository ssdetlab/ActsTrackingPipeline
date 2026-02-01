#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"

#include <cmath>
#include <memory>

#include "TrackingPipeline/Geometry/detail/BinningValueUtils.hpp"
#include "TrackingPipeline/Geometry/detail/SurfaceParameters.hpp"

/// TODO: cross-reference parameters with the paper

namespace E320Geometry {

using namespace Acts::UnitLiterals;

struct GeometryOptions {
  GeometryOptions() = default;
  ~GeometryOptions() = default;

  struct DipoleParameters {
    using AxisQuadruple =
        std::tuple<Acts::BinningValue, double, double, double>;

    DipoleParameters(AxisQuadruple tPrimary, AxisQuadruple tLong,
                     AxisQuadruple tShort)
        : rotationAnglePrimary(std::get<2>(tPrimary)),
          rotationAngleLong(std::get<2>(tLong)),
          rotationAngleShort(std::get<2>(tShort)) {
      center[detail::binningValueToIndex(std::get<0>(tPrimary))] =
          std::get<1>(tPrimary);
      center[detail::binningValueToIndex(std::get<0>(tLong))] =
          std::get<1>(tLong);
      center[detail::binningValueToIndex(std::get<0>(tShort))] =
          std::get<1>(tShort);

      field[detail::binningValueToIndex(std::get<0>(tPrimary))] =
          std::get<3>(tPrimary);
      field[detail::binningValueToIndex(std::get<0>(tLong))] =
          std::get<3>(tLong);
      field[detail::binningValueToIndex(std::get<0>(tShort))] =
          std::get<3>(tShort);
    }
    DipoleParameters() = delete;

    Acts::Vector3 center;

    double rotationAnglePrimary;
    double rotationAngleLong;
    double rotationAngleShort;

    Acts::Vector3 field;
  };

  using TrackingChamberParameters = std::vector<SurfaceParameters>;

  /// --------------------------------------------------------------
  /// General parameters

  /// Detector binning directions
  const Acts::BinningValue primaryBinValue = Acts::BinningValue::binX;
  const Acts::BinningValue longBinValue = Acts::BinningValue::binY;
  const Acts::BinningValue shortBinValue = Acts::BinningValue::binZ;

  const std::size_t primaryIdx = detail::binningValueToIndex(primaryBinValue);
  const std::size_t longIdx = detail::binningValueToIndex(longBinValue);
  const std::size_t shortIdx = detail::binningValueToIndex(shortBinValue);

  const Acts::Vector3 primaryDir =
      detail::binningValueToDirection(primaryBinValue);
  const Acts::Vector3 longDir = detail::binningValueToDirection(longBinValue);
  const Acts::Vector3 shortDir = detail::binningValueToDirection(shortBinValue);

  const std::size_t gapVolumeIdPrefactor = 50;
  const std::size_t magVolumeIdPrefactor = 100;
  const std::size_t ipVolumeIdPrefactor = 150;

  const double worldHalfLong = 2_m;
  const double worldHalfShort = 2_m;

  /// Rotation of the plane surfaces into the global frame
  const double toWorldAngleX = 0;
  const double toWorldAngleY = -M_PI_2;
  const double toWorldAngleZ = M_PI_2;

  /// --------------------------------------------------------------
  /// Parameters of the Be window

  const double beWindowHalfX = worldHalfLong;
  const double beWindowHalfY = worldHalfShort;

  const double beWindowThickness = 1_mm;

  /// --------------------------------------------------------------
  /// Parameters of the quads

  const double quad1HalfPrimary = 486.664_mm;
  const double quad1HalfLong = 32.33_mm;
  const double quad1HalfShort = 32.33_mm;
  const double quad1Gradient = -0.7637_T / 1_m;

  const double quad2HalfPrimary = 486.664_mm;
  const double quad2HalfLong = 32.33_mm;
  const double quad2HalfShort = 32.33_mm;
  const double quad2Gradient = 2.855_T / 1_m;

  const double quad3HalfPrimary = 486.664_mm;
  const double quad3HalfLong = 32.33_mm;
  const double quad3HalfShort = 32.33_mm;
  const double quad3Gradient = -0.7637_T / 1_m;

  /// --------------------------------------------------------------
  /// Parameters of the x-corrector

  const double xCorrectorHalfPrimary = 116.84_mm;
  const double xCorrectorHalfLong = 40_mm;
  const double xCorrectorHalfShort = 40_mm;

  const double xCorrectorFieldStrength = 0.026107_T;

  const double xCorrectorFieldPrimary = 0;
  const double xCorrectorFieldLong = xCorrectorFieldStrength;
  const double xCorrectorFieldShort = 0;

  /// --------------------------------------------------------------
  /// Parameters of the dipole

  const double dipoleHalfPrimary = 457.2_mm;
  const double dipoleHalfLong = 50.927_mm;
  const double dipoleHalfShort = 22.352_mm;

  const double dipoleFieldStrength = 0.2192_T;

  const double dipoleFieldPrimary = 0;
  const double dipoleFieldLong = 0;
  const double dipoleFieldShort = dipoleFieldStrength;

  /// --------------------------------------------------------------
  /// Parameters of the PDC window

  const double pdcWindowHalfX = 60_mm;
  const double pdcWindowHalfY = 25_mm;

  const double pdcWindowThickness = 0.51_mm;

  /// --------------------------------------------------------------
  /// Parameters of the tracking chambers

  /// Chip size in chips's local coordinates
  const double chipHalfX = 14.97088_mm;
  const double chipHalfY = 6.88128_mm;

  const double pixelHalfX = 14.62_um;
  const double pixelHalfY = 13.44_um;

  const double pixelThickness = 25_um;

  /// Volume spacing around the chips
  const double chipVolumeHalfSpacing = 1_mm;

  /// Distances within the detector box
  const double tcWindowToFirstChipDistance = 8.9_mm;
  const double tcWindowToLastChipDistance = 8.7_mm;
  const double interChipDistance = 20_mm;

  /// Transverse volume parameters
  const double tcHalfLong = chipHalfX;
  const double tcHalfShort = chipHalfY;

  /// --------------------------------------------------------------
  /// Be window placement

  // const double beWindowCenterPrimary = -842_mm;
  const double beWindowCenterPrimary = 0.01_um;
  const double beWindowCenterLong = 0_mm;
  const double beWindowCenterShort = 0_mm;

  const SurfaceParameters beWindowParameters{
      {primaryBinValue, beWindowCenterPrimary, toWorldAngleX},
      {longBinValue, beWindowCenterLong, toWorldAngleY},
      {shortBinValue, beWindowCenterShort, toWorldAngleZ},
      40};

  /// --------------------------------------------------------------
  /// BMP placement

  const double bpm1CenterPrimary = 4160_mm + 500_mm;
  const double bpm2CenterPrimary = 6390_mm + 500_mm;
  const double bpm3CenterPrimary = 8610_mm + 500_mm;

  const double bpmCenterLong = 0_mm;
  const double bpmCenterShort = 0_mm;

  const SurfaceParameters bpm1Parameters{
      {primaryBinValue, bpm1CenterPrimary, toWorldAngleX},
      {longBinValue, bpmCenterLong, toWorldAngleY},
      {shortBinValue, bpmCenterShort, toWorldAngleZ},
      41};
  const SurfaceParameters bpm2Parameters{
      {primaryBinValue, bpm2CenterPrimary, toWorldAngleX},
      {longBinValue, bpmCenterLong, toWorldAngleY},
      {shortBinValue, bpmCenterShort, toWorldAngleZ},
      42};
  const SurfaceParameters bpm3Parameters{
      {primaryBinValue, bpm3CenterPrimary, toWorldAngleX},
      {longBinValue, bpmCenterLong, toWorldAngleY},
      {shortBinValue, bpmCenterShort, toWorldAngleZ},
      43};

  /// --------------------------------------------------------------
  /// Quads placement

  const double quad1CenterPrimary = 4160_mm;
  const double quad1CenterLong = 0_mm;
  const double quad1CenterShort = 0_mm;

  const double quad2CenterPrimary = 6390_mm;
  const double quad2CenterLong = 0_mm;
  const double quad2CenterShort = 0_mm;

  const double quad3CenterPrimary = 8610_mm;
  const double quad3CenterLong = 0_mm;
  const double quad3CenterShort = 0_mm;

  /// --------------------------------------------------------------
  /// X-corrector placement

  const double xCorrectorCenterPrimary = 9994.63_mm;
  const double xCorrectorCenterLong = 0_mm;
  const double xCorrectorCenterShort = 0_mm;

  /// --------------------------------------------------------------
  /// Dipole placement

  const double dipoleCenterPrimary = 13060.61_mm;
  const double dipoleCenterLong = 0_mm;
  const double dipoleCenterShort = 0_mm;

  /// --------------------------------------------------------------
  /// PDC window placement

  const double pdcWindowCenterPrimary = 16549.745_mm;
  const double pdcWindowCenterLong = 118_mm;
  const double pdcWindowCenterShort = 0_mm;

  const SurfaceParameters pdcWindowParameters{
      {primaryBinValue, pdcWindowCenterPrimary, toWorldAngleX},
      {longBinValue, pdcWindowCenterLong, toWorldAngleY},
      {shortBinValue, pdcWindowCenterShort, toWorldAngleZ},
      30};

  /// --------------------------------------------------------------
  /// Tracking chamber placement

  const double ipTcDistance = 16674.4_mm;

  const double tcCenterLong = 93.775_mm;
  const double tcCenterShort = -0.61872_mm;

  const std::vector<SurfaceParameters> tcParameters{
      SurfaceParameters({primaryBinValue, ipTcDistance + 0 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tcCenterLong, toWorldAngleY},
                        {shortBinValue, tcCenterShort, toWorldAngleZ}, 10),
      SurfaceParameters{{primaryBinValue, ipTcDistance + 1 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tcCenterLong, toWorldAngleY},
                        {shortBinValue, tcCenterShort, toWorldAngleZ},
                        12},
      SurfaceParameters{{primaryBinValue, ipTcDistance + 2 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tcCenterLong, toWorldAngleY},
                        {shortBinValue, tcCenterShort, toWorldAngleZ},
                        14},
      SurfaceParameters{{primaryBinValue, ipTcDistance + 3 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tcCenterLong, toWorldAngleY},
                        {shortBinValue, tcCenterShort, toWorldAngleZ},
                        16},
      SurfaceParameters{{primaryBinValue, ipTcDistance + 4 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tcCenterLong, toWorldAngleY},
                        {shortBinValue, tcCenterShort, toWorldAngleZ},
                        18}};

  const double tcHalfPrimary =
      interChipDistance * (tcParameters.size() - 1) / 2.0;

  const double tcCenterPrimary = ipTcDistance + tcHalfPrimary;

  static const std::unique_ptr<const GeometryOptions>& instance() {
    if (!m_instance) {
      m_instance = std::make_unique<GeometryOptions>();
    }
    return m_instance;
  }

 protected:
  static std::unique_ptr<const GeometryOptions> m_instance;
};

}  // namespace E320Geometry
