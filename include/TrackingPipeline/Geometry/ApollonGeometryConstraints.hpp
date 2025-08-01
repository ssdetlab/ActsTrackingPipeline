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

  struct ChipParameters {
    using AxisTriple = std::tuple<Acts::BinningValue, double, double>;

    ChipParameters(AxisTriple tPrimary, AxisTriple tLong, AxisTriple tShort,
                   int id)
        : rotationAnglePrimary(std::get<2>(tPrimary)),
          rotationAngleLong(std::get<2>(tLong)),
          rotationAngleShort(std::get<2>(tShort)),
          geoId(id) {
      center[detail::binningValueToIndex(std::get<0>(tPrimary))] =
          std::get<1>(tPrimary);
      center[detail::binningValueToIndex(std::get<0>(tLong))] =
          std::get<1>(tLong);
      center[detail::binningValueToIndex(std::get<0>(tShort))] =
          std::get<1>(tShort);
    }
    ChipParameters() = delete;

    Acts::Vector3 center;

    double rotationAnglePrimary;
    double rotationAngleLong;
    double rotationAngleShort;

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

  using TrackingChamberParameters = std::vector<ChipParameters>;

  /// --------------------------------------------------------------
  /// General parameters

  /// Detector binning directions
  const Acts::BinningValue primaryBinValue = Acts::BinningValue::binX;
  const Acts::BinningValue longBinValue = Acts::BinningValue::binY;
  const Acts::BinningValue shortBinValue = Acts::BinningValue::binZ;

  const std::size_t gapVolumeIdPrefactor = 100;
  const std::size_t magVolumeIdPrefactor = 200;

  /// --------------------------------------------------------------
  /// Parameters of the tracking chambers

  /// Chip size in chips's local coordinates
  const double chipHalfX = 14.97088_mm;
  const double chipHalfY = 6.88128_mm;

  const double pixelHalfX = 14.62_um;
  const double pixelHalfY = 13.44_um;

  /// Volume spacing around the chips
  const double chipVolumeHalfSpacing = 1_mm;

  /// Rotation of the chips into the global frame
  const double toWorldAngleX = 0;
  const double toWorldAngleY = M_PI_2;
  const double toWorldAngleZ = M_PI_2;

  const double ipTc1Distance = 0_mm;
  const double interChipDistance = 10_mm;

  const double tcHalfLong = chipHalfX + chipVolumeHalfSpacing;
  const double tcHalfShort = chipHalfY + chipVolumeHalfSpacing;

  /// --------------------------------------------------------------
  /// Parameters of the dipole

  const double dipoleHalfPrimary = 60.5_mm;
  const double dipoleHalfLong = 30_mm;
  const double dipoleHalfShort = 10_mm;

  const double dipoleFieldPrimary = 0;
  const double dipoleFieldLong = 0;
  const double dipoleFieldShort = -0.35_T;

  /// --------------------------------------------------------------
  /// First tracking chamber placement

  const double tc1CenterLong = 0;
  const double tc1CenterShort = 0;

  const std::vector<ChipParameters> tc1Parameters{
      ChipParameters({primaryBinValue, ipTc1Distance + 0 * interChipDistance,
                      toWorldAngleX},
                     {longBinValue, tc1CenterLong, toWorldAngleY},
                     {shortBinValue, tc1CenterShort, toWorldAngleZ}, 10),
      ChipParameters{{primaryBinValue, ipTc1Distance + 1 * interChipDistance,
                      toWorldAngleX},
                     {longBinValue, tc1CenterLong, toWorldAngleY},
                     {shortBinValue, tc1CenterShort, toWorldAngleZ},
                     11},
      ChipParameters{{primaryBinValue, ipTc1Distance + 2 * interChipDistance,
                      toWorldAngleX},
                     {longBinValue, tc1CenterLong, toWorldAngleY},
                     {shortBinValue, tc1CenterShort, toWorldAngleZ},
                     12},
      ChipParameters{{primaryBinValue, ipTc1Distance + 3 * interChipDistance,
                      toWorldAngleX},
                     {longBinValue, tc1CenterLong, toWorldAngleY},
                     {shortBinValue, tc1CenterShort, toWorldAngleZ},
                     13},
      ChipParameters{{primaryBinValue, ipTc1Distance + 4 * interChipDistance,
                      toWorldAngleX},
                     {longBinValue, tc1CenterLong, toWorldAngleY},
                     {shortBinValue, tc1CenterShort, toWorldAngleZ},
                     14}};

  /// --------------------------------------------------------------
  /// Dipole placement

  const double tc1DipoleDistance = 20_mm;
  const double ipDipoleDistance =
      ipTc1Distance + interChipDistance * (tc1Parameters.size() - 1) +
      tc1DipoleDistance;

  const double dipoleCenterPrimary = ipDipoleDistance + dipoleHalfPrimary;
  const double dipoleCenterLong = 0;
  const double dipoleCenterShort = 0;

  const DipoleParameters dipoleParameters{
      {primaryBinValue, dipoleCenterPrimary, 0, dipoleFieldPrimary},
      {longBinValue, 0, 0, 0},
      {shortBinValue, 0, 0, -0.35_T}};

  /// --------------------------------------------------------------
  /// Second tracking chamber placement

  const double dipoleTc2Distance = 20_mm;
  const double ipTc2Distance =
      ipDipoleDistance + 2 * dipoleHalfPrimary + dipoleTc2Distance;

  const double tc2CenterLong = -15_mm;
  const double tc2CenterShort = 0;

  const std::vector<ChipParameters> tc2Parameters{
      ChipParameters{{primaryBinValue, ipTc2Distance + 0 * interChipDistance,
                      toWorldAngleX},
                     {longBinValue, tc2CenterLong, toWorldAngleY},
                     {shortBinValue, tc2CenterShort, toWorldAngleZ},
                     20},
      ChipParameters{{primaryBinValue, ipTc2Distance + 1 * interChipDistance,
                      toWorldAngleX},
                     {longBinValue, tc2CenterLong, toWorldAngleY},
                     {shortBinValue, tc2CenterShort, toWorldAngleZ},
                     21},
      ChipParameters{{primaryBinValue, ipTc2Distance + 2 * interChipDistance,
                      toWorldAngleX},
                     {longBinValue, tc2CenterLong, toWorldAngleY},
                     {shortBinValue, tc2CenterShort, toWorldAngleZ},
                     22},
      ChipParameters{{primaryBinValue, ipTc2Distance + 3 * interChipDistance,
                      toWorldAngleX},
                     {longBinValue, tc2CenterLong, toWorldAngleY},
                     {shortBinValue, tc2CenterShort, toWorldAngleZ},
                     23},
      ChipParameters{{primaryBinValue, ipTc2Distance + 4 * interChipDistance,
                      toWorldAngleX},
                     {longBinValue, tc2CenterLong, toWorldAngleY},
                     {shortBinValue, tc2CenterShort, toWorldAngleZ},
                     24}};

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
