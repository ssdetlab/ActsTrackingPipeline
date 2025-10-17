#include "TrackingPipeline/Geometry/ApollonGeometry.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/detail/CuboidalDetectorHelper.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <string>

#include <unistd.h>

#include "TrackingPipeline/Alignment/AlignableDetectorElement.hpp"
#include "TrackingPipeline/Geometry/ApollonGeometryConstraints.hpp"
#include "TrackingPipeline/MagneticField/ConstantBoundedField.hpp"

namespace ApollonGeometry {

using go = GeometryOptions;

std::shared_ptr<const Acts::Experimental::Detector> buildDetector(
    const Acts::GeometryContext& gctx, bool insertReferenceSurface) {
  const auto& goInst = *go::instance();

  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
      detectorVolumes;
  std::size_t nVolumes =
      (insertReferenceSurface)
          ? goInst.tc1Parameters.size() + goInst.tc2Parameters.size() + 1
          : goInst.tc1Parameters.size() + goInst.tc2Parameters.size() + 3;
  detectorVolumes.reserve(nVolumes);

  std::vector<std::shared_ptr<Acts::DetectorElementBase>> detectorElements;
  std::size_t nDetElements =
      (insertReferenceSurface)
          ? goInst.tc1Parameters.size() + goInst.tc2Parameters.size()
          : goInst.tc1Parameters.size() + goInst.tc2Parameters.size() + 1;
  detectorElements.reserve(nDetElements);

  auto chipBounds = std::make_shared<Acts::RectangleBounds>(goInst.chipHalfX,
                                                            goInst.chipHalfY);
  auto vcWindowBounds = std::make_shared<Acts::RectangleBounds>(
      goInst.vcWindowHalfX, goInst.vcWindowHalfY);

  std::size_t gapVolumeCounter = 0;
  std::size_t magVolumeCounter = 0;

  double negLongEdge = std::min(
      {goInst.tc1CenterLong - goInst.tcHalfLong,
       goInst.tc2CenterLong - goInst.tcHalfLong,
       goInst.dipoleCenterLong - goInst.dipoleHalfLong, -goInst.worldHalfLong});
  double posLongEdge = std::max(
      {goInst.tc1CenterLong + goInst.tcHalfLong,
       goInst.tc2CenterLong + goInst.tcHalfLong,
       goInst.dipoleCenterLong + goInst.dipoleHalfLong, goInst.worldHalfLong});

  double negShortEdge =
      std::min({goInst.tc1CenterShort - goInst.tcHalfShort,
                goInst.tc2CenterShort - goInst.tcHalfShort,
                goInst.dipoleCenterShort - goInst.dipoleHalfShort,
                -goInst.worldHalfShort});
  double posShortEdge =
      std::max({goInst.tc1CenterShort + goInst.tcHalfShort,
                goInst.tc2CenterShort + goInst.tcHalfShort,
                goInst.dipoleCenterShort + goInst.dipoleHalfShort,
                goInst.worldHalfShort});

  Acts::Material silicon = Acts::Material::fromMolarDensity(
      9.370_cm, 46.52_cm, 28.0855, 14, (2.329 / 28.0855) * 1_mol / 1_cm3);
  Acts::MaterialSlab siliconSlab(silicon, goInst.pixelThickness);
  auto siSurfMaterial =
      std::make_shared<Acts::HomogeneousSurfaceMaterial>(siliconSlab);

  Acts::Material alminium = Acts::Material::fromMolarDensity(
      8.897_cm, 25.81_cm, 26.9815385, 13, (2.699 / 26.9815385) * 1_mol / 1_cm3);
  Acts::MaterialSlab aluminiumSlab(alminium, goInst.vcWindowThickness);
  auto alSurfMaterial =
      std::make_shared<Acts::HomogeneousSurfaceMaterial>(aluminiumSlab);

  double longVolumeSize = (posLongEdge - negLongEdge) / 2.0;
  double shortVolumeSize = (posShortEdge - negShortEdge) / 2.0;

  double longVolumeCenter = (posLongEdge + negLongEdge) / 2.0;
  double shortVolumeCenter = (posShortEdge + negShortEdge) / 2.0;

  auto constructVolume =
      [&](double halfPrimary, double centerPrimary,
          const std::string& namePrefix, std::size_t id,
          const std::vector<std::shared_ptr<Acts::Surface>>& surfaces) {
        Acts::Transform3 transform = Acts::Transform3::Identity();

        std::cout << "-----------------------------------------------\n";
        std::cout << namePrefix + std::to_string(id) << "\n";

        std::array<double, 3> volBoundsArray;
        volBoundsArray.at(goInst.primaryIdx) = halfPrimary;
        volBoundsArray.at(goInst.longIdx) = longVolumeSize;
        volBoundsArray.at(goInst.shortIdx) = shortVolumeSize;
        auto volBounds =
            std::make_unique<Acts::CuboidVolumeBounds>(volBoundsArray);

        Acts::Transform3 volTransform = Acts::Transform3::Identity();

        Acts::Vector3 volTranslation(0, 0, 0);
        volTranslation[goInst.primaryIdx] = centerPrimary;
        volTranslation[goInst.longIdx] = longVolumeCenter;
        volTranslation[goInst.shortIdx] = shortVolumeCenter;

        volTransform.translate(volTranslation);
        auto vol = Acts::Experimental::DetectorVolumeFactory::construct(
            Acts::Experimental::defaultPortalAndSubPortalGenerator(), gctx,
            namePrefix + std::to_string(id), volTransform, std::move(volBounds),
            surfaces, {}, Acts::Experimental::tryNoVolumes(),
            Acts::Experimental::tryAllPortalsAndSurfaces());

        Acts::GeometryIdentifier volGeoId;
        volGeoId.setVolume(id);
        vol->assignGeometryId(volGeoId);

        std::cout << vol->extent(gctx) << "\n";

        int portId = 1;
        for (auto& port : vol->portalPtrs()) {
          Acts::GeometryIdentifier portGeoId;
          portGeoId.setVolume(id);
          portGeoId.setApproach(portId);
          port->surface().assignGeometryId(id);
          portId++;
        }

        detectorVolumes.push_back(vol);
      };

  auto constructSurface =
      [&](const go::SurfaceParameters& pars,
          const std::shared_ptr<Acts::RectangleBounds>& bounds) {
        Acts::Transform3 surfTransform = Acts::Transform3::Identity();

        Acts::RotationMatrix3 surfToWorldRotationPrimary =
            Acts::AngleAxis3(pars.toWorldAnglePrimary, goInst.primaryDir)
                .toRotationMatrix();
        Acts::RotationMatrix3 surfToWorldRotationLong =
            Acts::AngleAxis3(pars.toWorldAngleLong, goInst.longDir)
                .toRotationMatrix();
        Acts::RotationMatrix3 surfToWorldRotationShort =
            Acts::AngleAxis3(pars.toWorldAngleShort, goInst.shortDir)
                .toRotationMatrix();

        surfTransform.translate(pars.toWorldTranslation);

        surfTransform.rotate(surfToWorldRotationPrimary);
        surfTransform.rotate(surfToWorldRotationLong);
        surfTransform.rotate(surfToWorldRotationShort);

        return Acts::Surface::makeShared<Acts::PlaneSurface>(surfTransform,
                                                             bounds);
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
          auto surf = constructSurface(parameters, chipBounds);

          surf->assignSurfaceMaterial(siSurfMaterial);

          Acts::GeometryIdentifier surfGeoId;
          surfGeoId.setSensitive(parameters.geoId);
          surf->assignGeometryId(surfGeoId);

          detectorElements.push_back(std::make_shared<AlignableDetectorElement>(
              surf, surf->transform(gctx)));
          surf->assignDetectorElement(*detectorElements.back());

          constructVolume(goInst.interChipDistance / 2.0,
                          parameters.toWorldTranslation[goInst.primaryIdx],
                          "sensVol", parameters.geoId, {surf});
        }
      };

  auto constructDipoleVolume = [&](const go::DipoleParameters& pars) {
    int id = goInst.magVolumeIdPrefactor + magVolumeCounter;
    constructVolume(goInst.dipoleHalfPrimary, pars.center[goInst.primaryIdx],
                    "magVol", id, {});
    magVolumeCounter++;
  };

  auto constructGapVolume =
      [&](const std::shared_ptr<Acts::Experimental::DetectorVolume>& vol1,
          const std::shared_ptr<Acts::Experimental::DetectorVolume>& vol2) {
        double primaryCenter =
            (vol2->transform(gctx).translation()[goInst.primaryIdx] -
             vol2->volumeBounds().values().at(goInst.primaryIdx) +
             vol1->transform(gctx).translation()[goInst.primaryIdx] +
             vol1->volumeBounds().values().at(goInst.primaryIdx)) /
            2.0;
        double primaryHalf =
            (vol2->transform(gctx).translation()[goInst.primaryIdx] -
             vol2->volumeBounds().values().at(goInst.primaryIdx) -
             vol1->transform(gctx).translation()[goInst.primaryIdx] -
             vol1->volumeBounds().values().at(goInst.primaryIdx)) /
            2.0;
        int id = goInst.gapVolumeIdPrefactor + gapVolumeCounter;
        constructVolume(primaryHalf, primaryCenter, "gapVol", id, {});
        gapVolumeCounter++;
      };

  auto constructVCWindow = [&](const go::SurfaceParameters& pars) {
    auto surf = constructSurface(pars, vcWindowBounds);

    surf->assignSurfaceMaterial(alSurfMaterial);

    Acts::GeometryIdentifier surfGeoId;
    surfGeoId.setPassive(pars.geoId);
    surf->assignGeometryId(surfGeoId);

    constructVolume(goInst.chipVolumeHalfSpacing,
                    pars.toWorldTranslation[goInst.primaryIdx], "passVol",
                    pars.geoId, {surf});
  };

  if (insertReferenceSurface) {
    auto referenceSurfaceBounds = std::make_shared<Acts::RectangleBounds>(
        goInst.referenceSurfaceHalfX, goInst.referenceSurfaceHalfY);
    auto referenceSurface = constructSurface(goInst.referenceSurfaceParameters,
                                             referenceSurfaceBounds);
    Acts::GeometryIdentifier referenceSurfaceGeoId;
    referenceSurfaceGeoId.setSensitive(goInst.referenceSurfaceParameters.geoId);
    referenceSurface->assignGeometryId(referenceSurfaceGeoId);

    constructVolume(goInst.vcWindowCenterPrimary - goInst.chipVolumeHalfSpacing,
                    0, "ipVol", goInst.ipVolumeIdPrefactor, {referenceSurface});
  } else {
    constructVolume(goInst.vcWindowCenterPrimary - goInst.chipVolumeHalfSpacing,
                    0, "ipVol", goInst.ipVolumeIdPrefactor, {});
  }

  constructVCWindow(goInst.vcWindowParameters);
  std::size_t lastVCVolumeIdx = detectorVolumes.size() - 1;
  const auto lastVCVolume = detectorVolumes.at(lastVCVolumeIdx);

  constructTrackingChamber(goInst.tc1Parameters);
  std::size_t lastTc1VolumeIdx = detectorVolumes.size() - 1;
  const auto firstTc1Volume = detectorVolumes.at(lastVCVolumeIdx + 1);
  const auto lastTc1Volume = detectorVolumes.at(lastTc1VolumeIdx);

  constructDipoleVolume(goInst.dipoleParameters);
  std::size_t dipoleVolumeIdx = detectorVolumes.size() - 1;
  const auto dipoleVolume = detectorVolumes.at(dipoleVolumeIdx);

  constructTrackingChamber(goInst.tc2Parameters);
  const auto firstTc2Volume = detectorVolumes.at(dipoleVolumeIdx + 1);

  constructGapVolume(lastVCVolume, firstTc1Volume);
  constructGapVolume(lastTc1Volume, dipoleVolume);
  constructGapVolume(dipoleVolume, firstTc2Volume);

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

  return std::make_shared<ConstantBoundedField>(goInst.dipoleParameters.field,
                                                dipoleExtent);
}

}  // namespace ApollonGeometry
