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

#include "TrackingPipeline/Geometry/ApollonGeometryConstraints.hpp"
#include "TrackingPipeline/MagneticField/ConstantBoundedField.hpp"

namespace ApollonGeometry {

using go = GeometryOptions;

std::shared_ptr<const Acts::Experimental::Detector> buildDetector(
    const Acts::GeometryContext& gctx) {
  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
      detectorVolumes;
  detectorVolumes.reserve(go::instance()->tc1Parameters.size() +
                          go::instance()->tc2Parameters.size() + 1);
  auto chipBounds = std::make_shared<Acts::RectangleBounds>(
      go::instance()->chipHalfX, go::instance()->chipHalfY);
  auto vcWindowBounds = std::make_shared<Acts::RectangleBounds>(
      go::instance()->vcWindowHalfX, go::instance()->vcWindowHalfY);

  std::size_t gapVolumeCounter = 0;
  std::size_t magVolumeCounter = 0;

  double negLongEdge = std::min(
      {go::instance()->tc1CenterLong - go::instance()->tcHalfLong,
       go::instance()->tc2CenterLong - go::instance()->tcHalfLong,
       go::instance()->dipoleCenterLong - go::instance()->dipoleHalfLong,
       -go::instance()->worldHalfLong});
  double posLongEdge = std::max(
      {go::instance()->tc1CenterLong + go::instance()->tcHalfLong,
       go::instance()->tc2CenterLong + go::instance()->tcHalfLong,
       go::instance()->dipoleCenterLong + go::instance()->dipoleHalfLong,
       go::instance()->worldHalfLong});

  double negShortEdge = std::min(
      {go::instance()->tc1CenterShort - go::instance()->tcHalfShort,
       go::instance()->tc2CenterShort - go::instance()->tcHalfShort,
       go::instance()->dipoleCenterShort - go::instance()->dipoleHalfShort,
       -go::instance()->worldHalfShort});
  double posShortEdge = std::max(
      {go::instance()->tc1CenterShort + go::instance()->tcHalfShort,
       go::instance()->tc2CenterShort + go::instance()->tcHalfShort,
       go::instance()->dipoleCenterShort + go::instance()->dipoleHalfShort,
       go::instance()->worldHalfShort});

  Acts::Material silicon = Acts::Material::fromMolarDensity(
      9.370_cm, 46.52_cm, 28.0855, 14, (2.329 / 28.0855) * 1_mol / 1_cm3);
  Acts::MaterialSlab siliconSlab(silicon, go::instance()->pixelThickness);
  auto siSurfMaterial =
      std::make_shared<Acts::HomogeneousSurfaceMaterial>(siliconSlab);

  Acts::Material alminium = Acts::Material::fromMolarDensity(
      8.897_cm, 25.81_cm, 26.9815385, 13, (2.699 / 26.9815385) * 1_mol / 1_cm3);
  Acts::MaterialSlab aluminiumSlab(alminium, go::instance()->vcWindowThickness);
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
        volBoundsArray.at(go::instance()->primaryIdx) = halfPrimary;
        volBoundsArray.at(go::instance()->longIdx) = longVolumeSize;
        volBoundsArray.at(go::instance()->shortIdx) = shortVolumeSize;
        auto volBounds =
            std::make_unique<Acts::CuboidVolumeBounds>(volBoundsArray);

        Acts::Transform3 volTransform = Acts::Transform3::Identity();

        Acts::Vector3 volTranslation(0, 0, 0);
        volTranslation[go::instance()->primaryIdx] = centerPrimary;
        volTranslation[go::instance()->longIdx] = longVolumeCenter;
        volTranslation[go::instance()->shortIdx] = shortVolumeCenter;

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
            Acts::AngleAxis3(pars.toWorldAnglePrimary,
                             go::instance()->primaryDir)
                .toRotationMatrix();
        Acts::RotationMatrix3 surfToWorldRotationLong =
            Acts::AngleAxis3(pars.toWorldAngleLong, go::instance()->longDir)
                .toRotationMatrix();
        Acts::RotationMatrix3 surfToWorldRotationShort =
            Acts::AngleAxis3(pars.toWorldAngleShort, go::instance()->shortDir)
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
            (pars.front().toWorldTranslation[go::instance()->primaryIdx] -
             go::instance()->tcWindowToFirstChipDistance +
             pars.back().toWorldTranslation[go::instance()->primaryIdx] +
             go::instance()->tcWindowToLastChipDistance) /
            2.0;
        for (const auto& parameters : pars) {
          auto surf = constructSurface(parameters, chipBounds);

          surf->assignSurfaceMaterial(siSurfMaterial);

          Acts::GeometryIdentifier surfGeoId;
          surfGeoId.setSensitive(parameters.geoId);
          surf->assignGeometryId(surfGeoId);

          constructVolume(
              go::instance()->interChipDistance / 2.0,
              parameters.toWorldTranslation[go::instance()->primaryIdx],
              "sensVol", parameters.geoId, {surf});
        }
      };

  auto constructDipoleVolume = [&](const go::DipoleParameters& pars) {
    int id = go::instance()->magVolumeIdPrefactor + magVolumeCounter;
    constructVolume(go::instance()->dipoleHalfPrimary,
                    pars.center[go::instance()->primaryIdx], "magVol", id, {});
    magVolumeCounter++;
  };

  auto constructGapVolume =
      [&](const std::shared_ptr<Acts::Experimental::DetectorVolume>& vol1,
          const std::shared_ptr<Acts::Experimental::DetectorVolume>& vol2) {
        double primaryCenter =
            (vol2->transform(gctx).translation()[go::instance()->primaryIdx] -
             vol2->volumeBounds().values().at(go::instance()->primaryIdx) +
             vol1->transform(gctx).translation()[go::instance()->primaryIdx] +
             vol1->volumeBounds().values().at(go::instance()->primaryIdx)) /
            2.0;
        double primaryHalf =
            (vol2->transform(gctx).translation()[go::instance()->primaryIdx] -
             vol2->volumeBounds().values().at(go::instance()->primaryIdx) -
             vol1->transform(gctx).translation()[go::instance()->primaryIdx] -
             vol1->volumeBounds().values().at(go::instance()->primaryIdx)) /
            2.0;
        int id = go::instance()->gapVolumeIdPrefactor + gapVolumeCounter;
        constructVolume(primaryHalf, primaryCenter, "gapVol", id, {});
        gapVolumeCounter++;
      };

  auto constructVCWindow = [&](const go::SurfaceParameters& pars) {
    auto surf = constructSurface(pars, vcWindowBounds);

    surf->assignSurfaceMaterial(alSurfMaterial);

    Acts::GeometryIdentifier surfGeoId;
    surfGeoId.setPassive(pars.geoId);
    surf->assignGeometryId(surfGeoId);

    constructVolume(go::instance()->chipVolumeHalfSpacing,
                    pars.toWorldTranslation[go::instance()->primaryIdx],
                    "passVol", pars.geoId, {surf});
  };

  constructVolume(go::instance()->vcWindowCenterPrimary -
                      go::instance()->chipVolumeHalfSpacing,
                  0, "ipVol", go::instance()->ipVolumeIdPrefactor, {});

  constructVCWindow(go::instance()->vcWindowParameters);
  std::size_t lastVCVolumeIdx = detectorVolumes.size() - 1;
  const auto lastVCVolume = detectorVolumes.at(lastVCVolumeIdx);

  constructTrackingChamber(go::instance()->tc1Parameters);
  std::size_t lastTc1VolumeIdx = detectorVolumes.size() - 1;
  const auto firstTc1Volume = detectorVolumes.at(lastVCVolumeIdx + 1);
  const auto lastTc1Volume = detectorVolumes.at(lastTc1VolumeIdx);

  constructDipoleVolume(go::instance()->dipoleParameters);
  std::size_t dipoleVolumeIdx = detectorVolumes.size() - 1;
  const auto dipoleVolume = detectorVolumes.at(dipoleVolumeIdx);

  constructTrackingChamber(go::instance()->tc2Parameters);
  const auto firstTc2Volume = detectorVolumes.at(dipoleVolumeIdx + 1);

  constructGapVolume(lastVCVolume, firstTc1Volume);
  constructGapVolume(lastTc1Volume, dipoleVolume);
  constructGapVolume(dipoleVolume, firstTc2Volume);

  std::sort(
      detectorVolumes.begin(), detectorVolumes.end(),
      [&](const auto& vol1, const auto& vol2) {
        return (
            vol1->transform(gctx).translation()[go::instance()->primaryIdx] <
            vol2->transform(gctx).translation()[go::instance()->primaryIdx]);
      });

  auto portalContainer =
      Acts::Experimental::detail::CuboidalDetectorHelper::connect(
          gctx, detectorVolumes, go::instance()->primaryBinValue, {},
          Acts::Logging::VERBOSE);

  auto detector = Acts::Experimental::Detector::makeShared(
      "TelescopeDetector", detectorVolumes,
      Acts::Experimental::tryRootVolumes());

  return detector;
}

std::shared_ptr<Acts::MagneticFieldProvider> buildMagField(
    const Acts::GeometryContext& gctx) {
  Acts::Extent dipoleExtent;
  dipoleExtent.set(
      go::instance()->primaryBinValue,
      go::instance()->dipoleCenterPrimary - go::instance()->dipoleHalfPrimary,
      go::instance()->dipoleCenterPrimary + go::instance()->dipoleHalfPrimary);
  dipoleExtent.set(
      go::instance()->longBinValue,
      go::instance()->dipoleCenterLong - go::instance()->dipoleHalfLong,
      go::instance()->dipoleCenterLong + go::instance()->dipoleHalfLong);
  dipoleExtent.set(
      go::instance()->shortBinValue,
      go::instance()->dipoleCenterShort - go::instance()->dipoleHalfShort,
      go::instance()->dipoleCenterShort + go::instance()->dipoleHalfShort);

  return std::make_shared<ConstantBoundedField>(
      go::instance()->dipoleParameters.field, dipoleExtent);
}

}  // namespace ApollonGeometry
