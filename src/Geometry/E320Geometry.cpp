#include "TrackingPipeline/Geometry/E320Geometry.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/detail/CuboidalDetectorHelper.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <unistd.h>

#include "TrackingPipeline/Alignment/AlignableDetectorElement.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/detail/GeometryConstructionUtils.hpp"
#include "TrackingPipeline/MagneticField/CompositeMagField.hpp"
#include "TrackingPipeline/MagneticField/ConstantBoundedField.hpp"
#include "TrackingPipeline/MagneticField/IdealQuadrupoleMagField.hpp"

namespace E320Geometry {

using go = GeometryOptions;

std::shared_ptr<const Acts::Experimental::Detector> buildDetector(
    const Acts::GeometryContext& gctx) {
  const auto& goInst = *go::instance();

  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
      detectorVolumes;

  std::vector<std::shared_ptr<Acts::DetectorElementBase>> detectorElements;
  std::size_t nDetElements = goInst.tcParameters.size();
  detectorElements.reserve(nDetElements);

  auto chipBounds = std::make_shared<Acts::RectangleBounds>(goInst.chipHalfX,
                                                            goInst.chipHalfY);
  auto pdcWindowBounds = std::make_shared<Acts::RectangleBounds>(
      goInst.pdcWindowHalfX, goInst.pdcWindowHalfY);

  std::size_t gapVolumeCounter = 0;
  std::size_t magVolumeCounter = 0;

  double negLongEdge = std::min(
      {goInst.tcCenterLong - goInst.tcHalfLong,
       goInst.dipoleCenterLong - goInst.dipoleHalfLong, -goInst.worldHalfLong});
  double posLongEdge = std::max(
      {goInst.tcCenterLong + goInst.tcHalfLong,
       goInst.dipoleCenterLong + goInst.dipoleHalfLong, goInst.worldHalfLong});

  double negShortEdge =
      std::min({goInst.tcCenterShort - goInst.tcHalfShort,
                goInst.dipoleCenterShort - goInst.dipoleHalfShort,
                -goInst.worldHalfShort});
  double posShortEdge =
      std::max({goInst.tcCenterShort + goInst.tcHalfShort,
                goInst.dipoleCenterShort + goInst.dipoleHalfShort,
                goInst.worldHalfShort});

  Acts::Material silicon = Acts::Material::fromMolarDensity(
      9.370_cm, 46.52_cm, 28.0855, 14, (2.329 / 28.0855) * 1_mol / 1_cm3);
  Acts::MaterialSlab siliconSlab(silicon, goInst.pixelThickness);
  auto siSurfMaterial =
      std::make_shared<Acts::HomogeneousSurfaceMaterial>(siliconSlab);

  Acts::Material alminium = Acts::Material::fromMolarDensity(
      8.897_cm, 25.81_cm, 26.9815385, 13, (2.699 / 26.9815385) * 1_mol / 1_cm3);
  Acts::MaterialSlab aluminiumSlab(alminium, goInst.pdcWindowThickness);
  auto alSurfMaterial =
      std::make_shared<Acts::HomogeneousSurfaceMaterial>(aluminiumSlab);

  double longVolumeSize = (posLongEdge - negLongEdge) / 2.0;
  double shortVolumeSize = (posShortEdge - negShortEdge) / 2.0;

  double longVolumeCenter = (posLongEdge + negLongEdge) / 2.0;
  double shortVolumeCenter = (posShortEdge + negShortEdge) / 2.0;

  auto constructE320Volume =
      [&](double halfPrimary, double centerPrimary,
          const std::string& namePrefix, std::size_t id,
          const std::vector<std::shared_ptr<Acts::Surface>>& surfaces) {
        detectorVolumes.push_back(constructVolume(
            halfPrimary, longVolumeSize, shortVolumeSize, centerPrimary,
            longVolumeCenter, shortVolumeCenter, goInst.primaryIdx,
            goInst.longIdx, goInst.shortIdx, namePrefix, id, surfaces, gctx));
      };

  auto constructE320Surface =
      [&](const SurfaceParameters& pars,
          const std::shared_ptr<Acts::RectangleBounds>& bounds) {
        return constructSurface(pars, bounds, goInst.primaryDir, goInst.longDir,
                                goInst.shortDir);
      };

  auto constructTrackingChamber =
      [&](const go::TrackingChamberParameters& pars) {
        double tcBoxCenterPrimary =
            (pars.front().toWorldTranslation[goInst.primaryIdx] -
             goInst.tcWindowToFirstChipDistance +
             pars.back().toWorldTranslation[goInst.primaryIdx] +
             goInst.tcWindowToLastChipDistance) /
            2.0;
        for (const auto& parameters : pars) {
          auto surf = constructE320Surface(parameters, chipBounds);

          surf->assignSurfaceMaterial(siSurfMaterial);

          Acts::GeometryIdentifier surfGeoId;
          surfGeoId.setSensitive(parameters.geoId);
          surf->assignGeometryId(surfGeoId);

          detectorElements.push_back(std::make_shared<AlignableDetectorElement>(
              surf, surf->transform(gctx)));
          surf->assignDetectorElement(*detectorElements.back());

          constructE320Volume(goInst.interChipDistance / 2.0,
                              parameters.toWorldTranslation[goInst.primaryIdx],
                              "sensVol", parameters.geoId, {surf});
        }
      };

  auto constructMagVolume =
      [&](double halfPrimary, double centerPrimary,
          const std::string& namePrefix,
          const std::vector<std::shared_ptr<Acts::Surface>>& surfaces) {
        int id = goInst.magVolumeIdPrefactor + magVolumeCounter;
        constructE320Volume(halfPrimary, centerPrimary, namePrefix, id,
                            surfaces);
        magVolumeCounter++;
      };

  auto constructE320GapVolume =
      [&](const std::shared_ptr<Acts::Experimental::DetectorVolume>& vol1,
          const std::shared_ptr<Acts::Experimental::DetectorVolume>& vol2) {
        std::size_t id = goInst.gapVolumeIdPrefactor + gapVolumeCounter;
        detectorVolumes.push_back(
            constructGapVolume(vol1, vol2, goInst.primaryIdx, goInst.longIdx,
                               goInst.shortIdx, id, gctx));
        gapVolumeCounter++;
      };

  auto constructPDCWindow = [&](const SurfaceParameters& pars) {
    auto surf = constructE320Surface(pars, pdcWindowBounds);

    surf->assignSurfaceMaterial(alSurfMaterial);

    Acts::GeometryIdentifier surfGeoId;
    surfGeoId.setPassive(pars.geoId);
    surf->assignGeometryId(surfGeoId);

    constructE320Volume(goInst.chipVolumeHalfSpacing,
                        pars.toWorldTranslation[goInst.primaryIdx], "pdcVol",
                        pars.geoId, {surf});
  };

  // IP
  constructE320Volume(goInst.quad1CenterPrimary - goInst.quad1HalfPrimary, 0,
                      "ipVol", goInst.ipVolumeIdPrefactor, {});

  // Quad 1
  constructMagVolume(goInst.quad1HalfPrimary, goInst.quad1CenterPrimary,
                     "quad1", {});
  const auto& quad1Volume = detectorVolumes.back();

  // Quad 2
  constructMagVolume(goInst.quad2HalfPrimary, goInst.quad2CenterPrimary,
                     "quad2", {});
  const auto& quad2Volume = detectorVolumes.back();
  constructE320GapVolume(quad1Volume, quad2Volume);

  // Quad 3
  constructMagVolume(goInst.quad3HalfPrimary, goInst.quad3CenterPrimary,
                     "quad3", {});
  const auto quad3Volume = detectorVolumes.back();
  constructE320GapVolume(quad2Volume, quad3Volume);

  // X-Corrector
  constructMagVolume(goInst.xCorrectorHalfPrimary,
                     goInst.xCorrectorCenterPrimary, "xCorrector", {});
  const auto& xCorrectorVolume = detectorVolumes.back();
  constructE320GapVolume(quad3Volume, xCorrectorVolume);

  // Dipole
  constructMagVolume(goInst.dipoleHalfPrimary, goInst.dipoleCenterPrimary,
                     "dipole", {});
  const auto& dipoleVolume = detectorVolumes.back();
  constructE320GapVolume(xCorrectorVolume, dipoleVolume);

  // PDC window
  constructPDCWindow(goInst.pdcWindowParameters);
  std::size_t pdcVolumeIdx = detectorVolumes.size() - 1;
  const auto& pdcVolume = detectorVolumes.at(pdcVolumeIdx);
  constructE320GapVolume(dipoleVolume, pdcVolume);

  // Tracking chamber
  constructTrackingChamber(goInst.tcParameters);
  const auto& firstTcVolume = detectorVolumes.at(pdcVolumeIdx + 2);
  constructE320GapVolume(pdcVolume, firstTcVolume);

  std::sort(detectorVolumes.begin(), detectorVolumes.end(),
            [&](const auto& vol1, const auto& vol2) {
              return (vol1->transform(gctx).translation()[goInst.primaryIdx] <
                      vol2->transform(gctx).translation()[goInst.primaryIdx]);
            });

  auto portalContainer =
      Acts::Experimental::detail::CuboidalDetectorHelper::connect(
          gctx, detectorVolumes, goInst.primaryBinValue, {},
          Acts::Logging::VERBOSE);

  auto detector = Acts::Experimental::Detector::makeShared(
      "TelescopeDetector", detectorVolumes,
      Acts::Experimental::tryRootVolumes(), detectorElements);

  return detector;
}

std::shared_ptr<Acts::MagneticFieldProvider> buildMagField(
    const Acts::GeometryContext& gctx) {
  const auto& goInst = *go::instance();

  Acts::Extent quad1Extent;
  quad1Extent.set(goInst.primaryBinValue,
                  goInst.quad1CenterPrimary - goInst.quad1HalfPrimary,
                  goInst.quad1CenterPrimary + goInst.quad1HalfPrimary);
  quad1Extent.set(goInst.longBinValue,
                  goInst.quad1CenterLong - goInst.quad1HalfLong,
                  goInst.quad1CenterLong + goInst.quad1HalfLong);
  quad1Extent.set(goInst.shortBinValue,
                  goInst.quad1CenterShort - goInst.quad1HalfShort,
                  goInst.quad1CenterShort + goInst.quad1HalfShort);

  Acts::Extent quad2Extent;
  quad2Extent.set(goInst.primaryBinValue,
                  goInst.quad2CenterPrimary - goInst.quad2HalfPrimary,
                  goInst.quad2CenterPrimary + goInst.quad2HalfPrimary);
  quad2Extent.set(goInst.longBinValue,
                  goInst.quad2CenterLong - goInst.quad2HalfLong,
                  goInst.quad2CenterLong + goInst.quad2HalfLong);
  quad2Extent.set(goInst.shortBinValue,
                  goInst.quad2CenterShort - goInst.quad2HalfShort,
                  goInst.quad2CenterShort + goInst.quad2HalfShort);

  Acts::Extent quad3Extent;
  quad3Extent.set(goInst.primaryBinValue,
                  goInst.quad3CenterPrimary - goInst.quad3HalfPrimary,
                  goInst.quad3CenterPrimary + goInst.quad3HalfPrimary);
  quad3Extent.set(goInst.longBinValue,
                  goInst.quad3CenterLong - goInst.quad3HalfLong,
                  goInst.quad3CenterLong + goInst.quad3HalfLong);
  quad3Extent.set(goInst.shortBinValue,
                  goInst.quad3CenterShort - goInst.quad3HalfShort,
                  goInst.quad3CenterShort + goInst.quad3HalfShort);

  Acts::Extent xCorrectorExtent;
  xCorrectorExtent.set(
      goInst.primaryBinValue,
      goInst.xCorrectorCenterPrimary - goInst.xCorrectorHalfPrimary,
      goInst.xCorrectorCenterPrimary + goInst.xCorrectorHalfPrimary);
  xCorrectorExtent.set(goInst.longBinValue,
                       goInst.xCorrectorCenterLong - goInst.xCorrectorHalfLong,
                       goInst.xCorrectorCenterLong + goInst.xCorrectorHalfLong);
  xCorrectorExtent.set(
      goInst.shortBinValue,
      goInst.xCorrectorCenterShort - goInst.xCorrectorHalfShort,
      goInst.xCorrectorCenterShort + goInst.xCorrectorHalfShort);

  Acts::Extent dipoleExtent;
  dipoleExtent.set(goInst.primaryBinValue,
                   goInst.dipoleCenterPrimary - goInst.dipoleHalfPrimary,
                   goInst.dipoleCenterPrimary + goInst.dipoleHalfPrimary);
  dipoleExtent.set(goInst.longBinValue,
                   goInst.dipoleCenterLong - goInst.dipoleHalfLong,
                   goInst.dipoleCenterLong + goInst.dipoleHalfLong);
  dipoleExtent.set(goInst.shortBinValue,
                   goInst.dipoleCenterShort - goInst.dipoleHalfShort,
                   goInst.dipoleCenterShort + goInst.dipoleHalfShort);

  Acts::RotationMatrix3 quadRotation =
      Acts::AngleAxis3(-M_PI_2, Acts::Vector3::UnitY()).toRotationMatrix();

  Acts::Vector3 quad1Center;
  quad1Center[goInst.primaryIdx] = goInst.quad1CenterPrimary;
  quad1Center[goInst.longIdx] = goInst.quad1CenterLong;
  quad1Center[goInst.shortIdx] = goInst.quad1CenterShort;
  auto quad1Field = std::make_shared<IdealQuadrupoleMagField>(
      goInst.quad1Gradient, quad1Center, quadRotation);

  Acts::Vector3 quad2Center;
  quad2Center[goInst.primaryIdx] = goInst.quad2CenterPrimary;
  quad2Center[goInst.longIdx] = goInst.quad2CenterLong;
  quad2Center[goInst.shortIdx] = goInst.quad2CenterShort;
  auto quad2Field = std::make_shared<IdealQuadrupoleMagField>(
      goInst.quad2Gradient, quad2Center, quadRotation);

  Acts::Vector3 quad3Center;
  quad3Center[goInst.primaryIdx] = goInst.quad3CenterPrimary;
  quad3Center[goInst.longIdx] = goInst.quad3CenterLong;
  quad3Center[goInst.shortIdx] = goInst.quad3CenterShort;
  auto quad3Field = std::make_shared<IdealQuadrupoleMagField>(
      goInst.quad3Gradient, quad3Center, quadRotation);

  Acts::Vector3 xCorrectorB;
  xCorrectorB[goInst.primaryIdx] = goInst.xCorrectorFieldPrimary;
  xCorrectorB[goInst.longIdx] = goInst.xCorrectorFieldLong;
  xCorrectorB[goInst.shortIdx] = goInst.xCorrectorFieldShort;
  auto xCorrectorField =
      std::make_shared<ConstantBoundedField>(xCorrectorB, xCorrectorExtent);

  Acts::Vector3 dipoleB;
  dipoleB[goInst.primaryIdx] = goInst.dipoleFieldPrimary;
  dipoleB[goInst.longIdx] = goInst.dipoleFieldLong;
  dipoleB[goInst.shortIdx] = goInst.dipoleFieldShort;
  auto dipoleField =
      std::make_shared<ConstantBoundedField>(dipoleB, dipoleExtent);

  CompositeMagField::FieldComponents fieldComponents = {
      {quad1Extent, quad1Field},
      {quad2Extent, quad2Field},
      {quad3Extent, quad3Field},
      {xCorrectorExtent, xCorrectorField},
      {dipoleExtent, dipoleField}};

  return std::make_shared<CompositeMagField>(fieldComponents);
}

}  // namespace E320Geometry
