#include "TrackingPipeline/Geometry/ApollonGeometry.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/detail/CuboidalDetectorHelper.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
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

  std::size_t gapVolumeCounter = 0;
  std::size_t magVolumeCounter = 0;

  std::size_t primaryIdx =
      detail::binningValueToIndex(go::instance()->primaryBinValue);
  std::size_t longIdx =
      detail::binningValueToIndex(go::instance()->longBinValue);
  std::size_t shortIdx =
      detail::binningValueToIndex(go::instance()->shortBinValue);

  std::cout << "NEG LONG "
            << go::instance()->tc1CenterLong - go::instance()->tcHalfLong
            << " --- "
            << go::instance()->tc2CenterLong - go::instance()->tcHalfLong
            << " --- "
            << go::instance()->dipoleCenterLong - go::instance()->dipoleHalfLong
            << "\n";

  std::cout << "POS LONG "
            << go::instance()->tc1CenterLong + go::instance()->tcHalfLong
            << " --- "
            << go::instance()->tc2CenterLong + go::instance()->tcHalfLong
            << " --- "
            << go::instance()->dipoleCenterLong + go::instance()->dipoleHalfLong
            << "\n";
  double negLongEdge = std::min(
      {go::instance()->tc1CenterLong - go::instance()->tcHalfLong,
       go::instance()->tc2CenterLong - go::instance()->tcHalfLong,
       go::instance()->dipoleCenterLong - go::instance()->dipoleHalfLong});
  double posLongEdge = std::max(
      {go::instance()->tc1CenterLong + go::instance()->tcHalfLong,
       go::instance()->tc2CenterLong + go::instance()->tcHalfLong,
       go::instance()->dipoleCenterLong + go::instance()->dipoleHalfLong});

  double negShortEdge = std::min(
      {go::instance()->tc1CenterShort - go::instance()->tcHalfShort,
       go::instance()->tc2CenterShort - go::instance()->tcHalfShort,
       go::instance()->dipoleCenterShort - go::instance()->dipoleHalfShort});
  double posShortEdge = std::max(
      {go::instance()->tc1CenterShort + go::instance()->tcHalfShort,
       go::instance()->tc2CenterShort + go::instance()->tcHalfShort,
       go::instance()->dipoleCenterShort + go::instance()->dipoleHalfShort});

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
        volBoundsArray.at(primaryIdx) = halfPrimary;
        volBoundsArray.at(longIdx) = longVolumeSize;
        volBoundsArray.at(shortIdx) = shortVolumeSize;
        auto volBounds =
            std::make_unique<Acts::CuboidVolumeBounds>(volBoundsArray);

        Acts::Transform3 volTransform = Acts::Transform3::Identity();

        Acts::Vector3 volTranslation(0, 0, 0);
        volTranslation[primaryIdx] = centerPrimary;
        volTranslation[longIdx] = longVolumeCenter;
        volTranslation[shortIdx] = shortVolumeCenter;

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

  auto constructTrackingChamber =
      [&](const go::TrackingChamberParameters& pars) {
        double tcBoxCenterPrimary =
            (pars.front().toWorldTranslation[primaryIdx] -
             go::instance()->tcWindowToFirstChipDistance +
             pars.back().toWorldTranslation[primaryIdx] +
             go::instance()->tcWindowToLastChipDistance) /
            2.0;
        for (const auto& parameters : pars) {
          int id = parameters.geoId;

          Acts::Transform3 surfTransform = Acts::Transform3::Identity();

          Acts::RotationMatrix3 surfToWorldRotationPrimary =
              Acts::AngleAxis3(parameters.toWorldAnglePrimary,
                               detail::binningValueToDirection(
                                   go::instance()->primaryBinValue))
                  .toRotationMatrix();
          Acts::RotationMatrix3 surfToWorldRotationLong =
              Acts::AngleAxis3(
                  parameters.toWorldAngleLong,
                  detail::binningValueToDirection(go::instance()->longBinValue))
                  .toRotationMatrix();
          Acts::RotationMatrix3 surfToWorldRotationShort =
              Acts::AngleAxis3(parameters.toWorldAngleShort,
                               detail::binningValueToDirection(
                                   go::instance()->shortBinValue))
                  .toRotationMatrix();

          surfTransform.translate(parameters.toWorldTranslation);

          surfTransform.rotate(surfToWorldRotationPrimary);
          surfTransform.rotate(surfToWorldRotationLong);
          surfTransform.rotate(surfToWorldRotationShort);

          auto surf = Acts::Surface::makeShared<Acts::PlaneSurface>(
              surfTransform, chipBounds);

          Acts::GeometryIdentifier surfGeoId;
          surfGeoId.setSensitive(id);
          surf->assignGeometryId(surfGeoId);

          constructVolume(go::instance()->interChipDistance / 2.0,
                          parameters.toWorldTranslation[primaryIdx], "sensVol",
                          id, {surf});
        }
      };

  auto constructDipoleVolume = [&](const go::DipoleParameters& pars) {
    int id = go::instance()->magVolumeIdPrefactor + magVolumeCounter;
    constructVolume(go::instance()->dipoleHalfPrimary, pars.center[primaryIdx],
                    "magVol", id, {});
    magVolumeCounter++;
  };

  auto constructGapVolume =
      [&](const std::shared_ptr<Acts::Experimental::DetectorVolume>& vol1,
          const std::shared_ptr<Acts::Experimental::DetectorVolume>& vol2) {
        double primaryCenter =
            (vol2->transform(gctx).translation()[primaryIdx] -
             vol2->volumeBounds().values().at(primaryIdx) +
             vol1->transform(gctx).translation()[primaryIdx] +
             vol1->volumeBounds().values().at(primaryIdx)) /
            2.0;
        double primaryHalf = (vol2->transform(gctx).translation()[primaryIdx] -
                              vol2->volumeBounds().values().at(primaryIdx) -
                              vol1->transform(gctx).translation()[primaryIdx] -
                              vol1->volumeBounds().values().at(primaryIdx)) /
                             2.0;
        int id = go::instance()->gapVolumeIdPrefactor + gapVolumeCounter;
        constructVolume(primaryHalf, primaryCenter, "gapVol", id, {});
        gapVolumeCounter++;
      };

  if (go::instance()->ipTc1Distance > 0) {
    constructVolume(
        go::instance()->ipTc1Distance - go::instance()->interChipDistance / 2.0,
        0, "ipVol", go::instance()->ipVolumeIdPrefactor, {});
  }

  constructTrackingChamber(go::instance()->tc1Parameters);
  std::size_t lastTc1VolumeIdx = detectorVolumes.size() - 1;
  const auto lastTc1Volume = detectorVolumes.at(lastTc1VolumeIdx);

  constructDipoleVolume(go::instance()->dipoleParameters);
  std::size_t dipoleVolumeIdx = detectorVolumes.size() - 1;
  const auto dipoleVolume = detectorVolumes.at(dipoleVolumeIdx);

  constructTrackingChamber(go::instance()->tc2Parameters);
  const auto firstTc2Volume = detectorVolumes.at(dipoleVolumeIdx + 1);

  constructGapVolume(lastTc1Volume, dipoleVolume);
  constructGapVolume(dipoleVolume, firstTc2Volume);

  std::sort(detectorVolumes.begin(), detectorVolumes.end(),
            [&](const auto& vol1, const auto& vol2) {
              return (vol1->transform(gctx).translation()[primaryIdx] <
                      vol2->transform(gctx).translation()[primaryIdx]);
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
