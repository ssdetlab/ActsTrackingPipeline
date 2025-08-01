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

  std::cout << longVolumeCenter << "\n";
  std::cout << longVolumeSize << "\n";

  auto constructTrackingChamber =
      [&](const go::TrackingChamberParameters& pars) {
        for (const auto& parameters : pars) {
          int id = parameters.geoId;

          std::cout << "--------------------------------------------\n";

          Acts::Transform3 surfTransform = Acts::Transform3::Identity();

          Acts::RotationMatrix3 surfRotationPrimary =
              Acts::AngleAxis3(parameters.rotationAnglePrimary,
                               detail::binningValueToDirection(
                                   go::instance()->primaryBinValue))
                  .toRotationMatrix();
          Acts::RotationMatrix3 surfRotationLong =
              Acts::AngleAxis3(
                  parameters.rotationAngleLong,
                  detail::binningValueToDirection(go::instance()->longBinValue))
                  .toRotationMatrix();
          Acts::RotationMatrix3 surfRotationShort =
              Acts::AngleAxis3(parameters.rotationAngleShort,
                               detail::binningValueToDirection(
                                   go::instance()->shortBinValue))
                  .toRotationMatrix();

          surfTransform.translate(parameters.center);
          surfTransform.rotate(surfRotationPrimary);
          surfTransform.rotate(surfRotationLong);
          surfTransform.rotate(surfRotationShort);

          auto surf = Acts::Surface::makeShared<Acts::PlaneSurface>(
              surfTransform, chipBounds);

          Acts::GeometryIdentifier surfGeoId;
          surfGeoId.setSensitive(id);
          surf->assignGeometryId(surfGeoId);

          std::array<double, 3> volBoundsArray;
          volBoundsArray.at(primaryIdx) =
              go::instance()->interChipDistance / 2.0;
          volBoundsArray.at(longIdx) = longVolumeSize;
          volBoundsArray.at(shortIdx) = shortVolumeSize;
          auto volBounds =
              std::make_unique<Acts::CuboidVolumeBounds>(volBoundsArray);

          Acts::Transform3 volTransform = Acts::Transform3::Identity();

          Acts::Vector3 volTranslation(0, 0, 0);
          volTranslation[primaryIdx] = parameters.center[primaryIdx];
          volTranslation[longIdx] = longVolumeCenter;
          volTranslation[shortIdx] = shortVolumeCenter;

          volTransform.translate(volTranslation);
          auto vol = Acts::Experimental::DetectorVolumeFactory::construct(
              Acts::Experimental::defaultPortalAndSubPortalGenerator(), gctx,
              "sensVol" + std::to_string(id), volTransform,
              std::move(volBounds), {surf}, {},
              Acts::Experimental::tryNoVolumes(),
              Acts::Experimental::tryAllPortalsAndSurfaces());

          std::cout << vol->name() << "\n";
          std::cout << vol->extent(gctx) << "\n";

          Acts::GeometryIdentifier volGeoId;
          volGeoId.setVolume(id);
          vol->assignGeometryId(volGeoId);

          int portId = 1;
          for (auto& port : vol->portalPtrs()) {
            Acts::GeometryIdentifier portGeoId;
            portGeoId.setVolume(id);
            portGeoId.setApproach(portId);
            port->surface().assignGeometryId(id);
            portId++;
          }

          detectorVolumes.push_back(vol);
        }
      };

  auto constructDipoleVolume = [&](const go::DipoleParameters& pars) {
    Acts::Transform3 transform = Acts::Transform3::Identity();

    std::cout << "--------------------------------------------\n";

    std::array<double, 3> volBoundsArray;
    volBoundsArray.at(primaryIdx) = go::instance()->dipoleHalfPrimary;
    volBoundsArray.at(longIdx) = longVolumeSize;
    volBoundsArray.at(shortIdx) = shortVolumeSize;
    auto volBounds = std::make_unique<Acts::CuboidVolumeBounds>(volBoundsArray);

    int id = go::instance()->magVolumeIdPrefactor + magVolumeCounter;
    magVolumeCounter++;

    Acts::Transform3 volTransform = Acts::Transform3::Identity();

    Acts::Vector3 volTranslation(0, 0, 0);
    volTranslation[primaryIdx] = pars.center[primaryIdx];
    volTranslation[longIdx] = longVolumeCenter;
    volTranslation[shortIdx] = shortVolumeCenter;

    volTransform.translate(volTranslation);
    auto vol = Acts::Experimental::DetectorVolumeFactory::construct(
        Acts::Experimental::defaultPortalAndSubPortalGenerator(), gctx,
        "magVol" + std::to_string(id), volTransform, std::move(volBounds), {},
        {}, Acts::Experimental::tryNoVolumes(),
        Acts::Experimental::tryAllPortalsAndSurfaces());

    Acts::GeometryIdentifier volGeoId;
    volGeoId.setVolume(id);
    vol->assignGeometryId(volGeoId);

    std::cout << vol->name() << "\n";
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

  auto constructGapVolume =
      [&](const std::shared_ptr<Acts::Experimental::DetectorVolume>& vol1,
          const std::shared_ptr<Acts::Experimental::DetectorVolume>& vol2) {
        std::cout << "--------------------------------------------\n";

        Acts::Transform3 volTransform = Acts::Transform3::Identity();
        Acts::Vector3 translation(0, 0, 0);
        translation[primaryIdx] =
            (vol2->transform(gctx).translation()[primaryIdx] -
             vol2->volumeBounds().values().at(primaryIdx) +
             vol1->transform(gctx).translation()[primaryIdx] +
             vol1->volumeBounds().values().at(primaryIdx)) /
            2.0;
        translation[longIdx] = vol2->transform(gctx).translation()[longIdx];
        translation[shortIdx] = vol2->transform(gctx).translation()[shortIdx];
        volTransform.translate(translation);

        std::array<double, 3> volBoundsArray;
        volBoundsArray.at(primaryIdx) =
            (vol2->transform(gctx).translation()[primaryIdx] -
             vol2->volumeBounds().values().at(primaryIdx) -
             vol1->transform(gctx).translation()[primaryIdx] -
             vol1->volumeBounds().values().at(primaryIdx)) /
            2.0;
        volBoundsArray.at(longIdx) = vol1->volumeBounds().values().at(longIdx);
        volBoundsArray.at(shortIdx) =
            vol1->volumeBounds().values().at(shortIdx);
        auto volBounds =
            std::make_unique<Acts::CuboidVolumeBounds>(volBoundsArray);

        auto vol = Acts::Experimental::DetectorVolumeFactory::construct(
            Acts::Experimental::defaultPortalAndSubPortalGenerator(), gctx,
            "gapVol" + std::to_string(gapVolumeCounter), volTransform,
            std::move(volBounds), {}, {}, Acts::Experimental::tryNoVolumes(),
            Acts::Experimental::tryAllPortalsAndSurfaces());

        int id = go::instance()->gapVolumeIdPrefactor + gapVolumeCounter;
        gapVolumeCounter++;
        Acts::GeometryIdentifier volGeoId;
        volGeoId.setVolume(id);
        vol->assignGeometryId(volGeoId);

        int portId = 1;
        for (auto& port : vol->portalPtrs()) {
          Acts::GeometryIdentifier portGeoId;
          portGeoId.setVolume(id);
          portGeoId.setApproach(portId);
          port->surface().assignGeometryId(id);
          portId++;
        }

        std::cout << vol->name() << "\n";
        std::cout << vol->extent(gctx) << "\n";

        detectorVolumes.push_back(vol);
      };

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
