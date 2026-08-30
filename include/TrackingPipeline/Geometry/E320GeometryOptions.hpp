#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/BinUtility.hpp"

#include <cmath>
#include <cstddef>
#include <memory>

#include "TrackingPipeline/Geometry/detail/SurfaceParameters.hpp"

namespace E320 {

using namespace Acts::UnitLiterals;

struct GeometryOptions {
  GeometryOptions();
  ~GeometryOptions() = default;

  using TrackingChamberParameters = std::vector<SurfaceParameters>;

  /// --------------------------------------------------------------
  /// General parameters

  /// Detector binning directions
  Acts::BinningValue primaryBinValue = Acts::BinningValue::binX;
  Acts::BinningValue longBinValue = Acts::BinningValue::binY;
  Acts::BinningValue shortBinValue = Acts::BinningValue::binZ;

  std::size_t primaryIdx;
  std::size_t longIdx;
  std::size_t shortIdx;

  Acts::Vector3 primaryDir;
  Acts::Vector3 longDir;
  Acts::Vector3 shortDir;

  // Volume geometry ID prefactors
  std::size_t magVolumeIdPrefactor;
  std::size_t gapVolumeIdPrefactor;

  double worldHalfLong;
  double worldHalfShort;

  /// Rotation of the plane surfaces into the global frame
  double toWorldAngleX;
  double toWorldAngleY;
  double toWorldAngleZ;

  /// --------------------------------------------------------------
  /// Parameters of the BPMs

  double bpmHalfX;
  double bpmHalfY;

  double bpmThickness;

  double bpm0CenterPrimary;
  double bpm1CenterPrimary;
  double bpm2CenterPrimary;
  double bpm3CenterPrimary;

  double bpm0CenterLong;
  double bpm1CenterLong;
  double bpm2CenterLong;
  double bpm3CenterLong;

  double bpm0CenterShort;
  double bpm1CenterShort;
  double bpm2CenterShort;
  double bpm3CenterShort;

  SurfaceParameters bpm0Parameters;
  SurfaceParameters bpm1Parameters;
  SurfaceParameters bpm2Parameters;
  SurfaceParameters bpm3Parameters;

  /// --------------------------------------------------------------
  /// Parameters of the Be window

  double beWindowHalfX;
  double beWindowHalfY;

  double beWindowThickness;

  double beWindowCenterPrimary;
  double beWindowCenterLong;
  double beWindowCenterShort;

  SurfaceParameters beWindowParameters;

  /// --------------------------------------------------------------
  /// Parameters of the IP surface

  double ipSurfaceHalfX;
  double ipSurfaceHalfY;

  double ipSurfaceThickness;

  double ipSurfaceCenterPrimary;
  double ipSurfaceCenterLong;
  double ipSurfaceCenterShort;

  SurfaceParameters ipSurfaceParameters;

  /// --------------------------------------------------------------
  /// Parameters of the PDC window

  double pdcWindowHalfX;
  double pdcWindowHalfY;

  double pdcWindowThickness;

  /// Material binning
  Acts::BinUtility pdcWindowMaterialBinningX;
  Acts::BinUtility pdcWindowMaterialBinningY;

  double pdcWindowCenterPrimary;
  double pdcWindowCenterLong;
  double pdcWindowCenterShort;

  SurfaceParameters pdcWindowParameters;

  /// --------------------------------------------------------------
  /// Parameters of the composite field

  std::size_t compositeMagFieldId;

  /// --------------------------------------------------------------
  /// Parameters of the quads

  std::size_t quad1Id;
  std::size_t quad2Id;
  std::size_t quad3Id;

  double quad1HalfPrimary;
  double quad1HalfLong;
  double quad1HalfShort;

  double quad2HalfPrimary;
  double quad2HalfLong;
  double quad2HalfShort;

  double quad3HalfPrimary;
  double quad3HalfLong;
  double quad3HalfShort;

  double quad1Gradient;
  double quad2Gradient;
  double quad3Gradient;

  double quad0CenterPrimary;
  double quad0CenterLong;
  double quad0CenterShort;

  double quad1CenterPrimary;
  double quad1CenterLong;
  double quad1CenterShort;

  double quad2CenterPrimary;
  double quad2CenterLong;
  double quad2CenterShort;

  double quad3CenterPrimary;
  double quad3CenterLong;
  double quad3CenterShort;

  /// --------------------------------------------------------------
  /// Parameters of the x-corrector

  std::size_t xCorrectorId;
  double xCorrectorHalfPrimary;
  double xCorrectorHalfLong;
  double xCorrectorHalfShort;

  double xCorrectorFieldStrength;

  std::size_t xCorrectorDirIdx;

  double xCorrectorFieldPrimary;
  double xCorrectorFieldLong;
  double xCorrectorFieldShort;

  double xCorrectorCenterPrimary;
  double xCorrectorCenterLong;
  double xCorrectorCenterShort;

  /// --------------------------------------------------------------
  /// Parameters of the dipole

  std::size_t dipoleId;
  double dipoleHalfPrimary;
  double dipoleHalfLong;
  double dipoleHalfShort;

  double dipoleFieldStrength;

  std::size_t dipoleDirIdx;

  double dipoleFieldPrimary;
  double dipoleFieldLong;
  double dipoleFieldShort;

  double dipoleCenterPrimary;
  double dipoleCenterLong;
  double dipoleCenterShort;

  /// --------------------------------------------------------------
  /// Parameters of the tracking chambers

  /// Chip size in chips's local coordinates
  double chipHalfX;
  double chipHalfY;

  double pixelHalfX;
  double pixelHalfY;

  double pixelThickness;

  /// Allpix2 measurement resolutions
  std::unordered_map<std::size_t, std::pair<double, double>> allpixErrors;

  /// Volume spacing around the chips
  double chipVolumeHalfSpacing;

  /// Distance between the chips
  double interChipDistance;
  double tcWindowToFirstChipDistance;
  double tcWindowToLastChipDistance;

  /// Transverse volume parameters
  double tcHalfLong;
  double tcHalfShort;

  /// Material binning
  Acts::BinUtility chipMaterialBinningX;
  Acts::BinUtility chipMaterialBinningY;

  // Placement
  double ipTcDistance;

  double tcCenterLong;
  double tcCenterShort;

  std::vector<SurfaceParameters> tcParameters;

  double tcHalfPrimary;
  double tcCenterPrimary;

  /// @brief static call for creating an instance
  static const std::unique_ptr<const GeometryOptions>& instance() {
    if (!m_instance) {
      m_instance = std::make_unique<GeometryOptions>();
    }
    return m_instance;
  }

 protected:
  /// Static instance pointer
  static std::unique_ptr<const GeometryOptions> m_instance;
};

}  // namespace E320
