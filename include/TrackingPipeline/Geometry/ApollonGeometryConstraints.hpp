#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace detail {

inline Acts::BinningValue directionToBinningValue(const Acts::Vector3& dir) {
  if (dir == Acts::Vector3::UnitX()) {
    return Acts::BinningValue::binX;
  } else if (dir == Acts::Vector3::UnitY()) {
    return Acts::BinningValue::binY;
  } else if (dir == Acts::Vector3::UnitZ()) {
    return Acts::BinningValue::binZ;
  } else {
    throw std::runtime_error("Invalid vector to binning value conversion");
  }
}

inline Acts::Vector3 binningValueToDirection(
    const Acts::BinningValue& binningValue) {
  switch (binningValue) {
    case Acts::BinningValue::binX:
      return Acts::Vector3::UnitX();
    case Acts::BinningValue::binY:
      return Acts::Vector3::UnitY();
    case Acts::BinningValue::binZ:
      return Acts::Vector3::UnitZ();
    default:
      throw std::runtime_error("Invalid binning value to vector conversion");
  }
}

inline std::size_t binningValueToIndex(const Acts::BinningValue& binningValue) {
  switch (binningValue) {
    case Acts::BinningValue::binX:
      return 0;
    case Acts::BinningValue::binY:
      return 1;
    case Acts::BinningValue::binZ:
      return 2;
    default:
      throw std::runtime_error("Invalid binning value to index conversion");
  }
}

};  // namespace detail

namespace ApollonGeometry {

using namespace Acts::UnitLiterals;

struct GeometryOptions {
  GeometryOptions() = default;
  ~GeometryOptions() = default;

  struct SurfaceParameters {
    using AxisPars = std::tuple<Acts::BinningValue, double, double>;

    SurfaceParameters(AxisPars tPrimary, AxisPars tLong, AxisPars tShort,
                      int id)
        : toWorldAnglePrimary(std::get<2>(tPrimary)),
          toWorldAngleLong(std::get<2>(tLong)),
          toWorldAngleShort(std::get<2>(tShort)),
          geoId(id) {
      toWorldTranslation[detail::binningValueToIndex(std::get<0>(tPrimary))] =
          std::get<1>(tPrimary);
      toWorldTranslation[detail::binningValueToIndex(std::get<0>(tLong))] =
          std::get<1>(tLong);
      toWorldTranslation[detail::binningValueToIndex(std::get<0>(tShort))] =
          std::get<1>(tShort);
    }
    SurfaceParameters() = delete;

    Acts::Vector3 toWorldTranslation;

    double toWorldAnglePrimary;
    double toWorldAngleLong;
    double toWorldAngleShort;

    int geoId;
  };

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

  const double setupRotationAngle = -2_degree;

  /// Rotation of the plane surfaces into the global frame
  const double toWorldAngleX = 0;
  const double toWorldAngleY = -M_PI_2;
  const double toWorldAngleZ = M_PI_2;

  /// --------------------------------------------------------------
  /// Parameters of the VC window

  const double vcWindowCenterPrimary = 1158.6_mm;

  const double vcWindowHalfX = 95_mm;
  const double vcWindowHalfY = 15_mm;

  const double vcWindowThickness = 0.5_mm;

  const double vcWindowAnglePrimary = M_PI_4;
  const double vcWindowAngleShort = 0;
  const double vcWindowAngleLong = 0;

  const SurfaceParameters vcWindowParameters{
      {primaryBinValue, vcWindowCenterPrimary, toWorldAngleX},
      {longBinValue, tc1CenterLong, toWorldAngleY},
      {shortBinValue, tc1CenterShort, toWorldAngleZ + vcWindowAnglePrimary},
      30};

  const double vcRad = vcWindowCenterPrimary + 6_mm;

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
  const double tcWindowToFirstChipDistance = 11.81_mm;
  const double tcWindowToLastChipDistance = 13.59_mm;
  const double interChipDistance = 20_mm;

  /// Transverse volume parameters
  const double tcHalfLong = chipHalfX + chipVolumeHalfSpacing;
  const double tcHalfShort = chipHalfY + chipVolumeHalfSpacing;

  /// --------------------------------------------------------------
  /// Parameters of the dipole

  const double dipoleAlCoverThickness = 2_mm;

  const double dipoleHalfPrimary = 60.5_mm;
  const double dipoleHalfLong = 30_mm;
  const double dipoleHalfShort = 10_mm;

  const double dipoleFieldPrimary = 0;
  const double dipoleFieldLong = 0;
  const double dipoleFieldShort = -0.35_T;

  /// --------------------------------------------------------------
  /// First tracking chamber placement

  const double vcExitTc1Distance = 20_mm + tcWindowToFirstChipDistance;
  const double ipTc1Distance = vcRad + vcExitTc1Distance;

  const double tc1CenterLong = 0_mm;
  const double tc1CenterShort = 0_mm;
  // const double tc1CenterLong = -50_mm;
  // const double tc1CenterShort = -50_mm;

  const std::vector<SurfaceParameters> tc1Parameters{
      SurfaceParameters({primaryBinValue, ipTc1Distance + 0 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tc1CenterLong, toWorldAngleY},
                        {shortBinValue, tc1CenterShort, toWorldAngleZ}, 10),
      SurfaceParameters{{primaryBinValue, ipTc1Distance + 1 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tc1CenterLong, toWorldAngleY},
                        {shortBinValue, tc1CenterShort, toWorldAngleZ},
                        12},
      SurfaceParameters{{primaryBinValue, ipTc1Distance + 2 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tc1CenterLong, toWorldAngleY},
                        {shortBinValue, tc1CenterShort, toWorldAngleZ},
                        14},
      SurfaceParameters{{primaryBinValue, ipTc1Distance + 3 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tc1CenterLong, toWorldAngleY},
                        {shortBinValue, tc1CenterShort, toWorldAngleZ},
                        16},
      SurfaceParameters{{primaryBinValue, ipTc1Distance + 4 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tc1CenterLong, toWorldAngleY},
                        {shortBinValue, tc1CenterShort, toWorldAngleZ},
                        18}};

  const double tc1HalfPrimary =
      interChipDistance * (tc1Parameters.size() - 1) / 2.0 +
      chipVolumeHalfSpacing;

  const double tc1CenterPrimary = ipTc1Distance + tc1HalfPrimary;

  /// --------------------------------------------------------------
  /// Dipole placement

  const double tc1DipoleDistance = 20_mm;
  const double ipDipoleDistance =
      ipTc1Distance + interChipDistance * (tc1Parameters.size() - 1) +
      tcWindowToLastChipDistance + tc1DipoleDistance + dipoleAlCoverThickness;

  const double dipoleCenterPrimary = ipDipoleDistance + dipoleHalfPrimary;
  // const double dipoleCenterLong = -50_mm;
  // const double dipoleCenterShort = -50_mm;
  const double dipoleCenterLong = 0_mm;
  const double dipoleCenterShort = 0_mm;

  const DipoleParameters dipoleParameters{
      {primaryBinValue, dipoleCenterPrimary, 0, dipoleFieldPrimary},
      {longBinValue, 0, 0, 0},
      {shortBinValue, 0, 0, -0.35_T}};

  /// --------------------------------------------------------------
  /// Second tracking chamber placement

  const double dipoleTc2Distance = 20_mm;
  const double ipTc2Distance = ipDipoleDistance + 2 * dipoleHalfPrimary +
                               dipoleAlCoverThickness + dipoleTc2Distance +
                               tcWindowToFirstChipDistance;

  // const double tc2CenterLong = -65_mm;
  // const double tc2CenterShort = -50_mm;
  const double tc2CenterLong = 0_mm;
  const double tc2CenterShort = 0_mm;

  const std::vector<SurfaceParameters> tc2Parameters{
      SurfaceParameters{{primaryBinValue, ipTc2Distance + 0 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tc2CenterLong, toWorldAngleY},
                        {shortBinValue, tc2CenterShort, toWorldAngleZ},
                        20},
      SurfaceParameters{{primaryBinValue, ipTc2Distance + 1 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tc2CenterLong, toWorldAngleY},
                        {shortBinValue, tc2CenterShort, toWorldAngleZ},
                        22},
      SurfaceParameters{{primaryBinValue, ipTc2Distance + 2 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tc2CenterLong, toWorldAngleY},
                        {shortBinValue, tc2CenterShort, toWorldAngleZ},
                        24},
      SurfaceParameters{{primaryBinValue, ipTc2Distance + 3 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tc2CenterLong, toWorldAngleY},
                        {shortBinValue, tc2CenterShort, toWorldAngleZ},
                        26},
      SurfaceParameters{{primaryBinValue, ipTc2Distance + 4 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tc2CenterLong, toWorldAngleY},
                        {shortBinValue, tc2CenterShort, toWorldAngleZ},
                        28}};

  const double tc2HalfPrimary =
      interChipDistance * (tc2Parameters.size() - 1) / 2.0 +
      chipVolumeHalfSpacing;

  const double tc2CenterPrimary = ipTc2Distance + tc2HalfPrimary;

  static const std::unique_ptr<const GeometryOptions>& instance() {
    if (!m_instance) {
      m_instance = std::make_unique<GeometryOptions>();
    }
    return m_instance;
  }

 protected:
  static std::unique_ptr<const GeometryOptions> m_instance;
};

}  // namespace ApollonGeometry
