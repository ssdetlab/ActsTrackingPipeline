#include "TrackingPipeline/Geometry/detail/GeometryConstructionUtils.hpp"

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>

std::shared_ptr<Acts::Experimental::DetectorVolume> constructVolume(
    double halfPrimary, double halfLong, double halfShort, double centerPrimary,
    double centerLong, double centerShort, std::size_t primaryIdx,
    std::size_t longIdx, std::size_t shortIdx, const std::string& namePrefix,
    std::size_t id, const std::vector<std::shared_ptr<Acts::Surface>>& surfaces,
    const Acts::GeometryContext& gctx) {
  Acts::Transform3 transform = Acts::Transform3::Identity();

  std::cout << "-----------------------------------------------\n";
  std::cout << namePrefix + std::to_string(id) << "\n";

  std::array<double, 3> volBoundsArray;
  volBoundsArray.at(primaryIdx) = halfPrimary;
  volBoundsArray.at(longIdx) = halfLong;
  volBoundsArray.at(shortIdx) = halfShort;
  auto volBounds = std::make_unique<Acts::CuboidVolumeBounds>(volBoundsArray);

  Acts::Transform3 volTransform = Acts::Transform3::Identity();

  Acts::Vector3 volTranslation(0, 0, 0);
  volTranslation[primaryIdx] = centerPrimary;
  volTranslation[longIdx] = centerLong;
  volTranslation[shortIdx] = centerShort;

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

  return vol;
};

std::shared_ptr<Acts::Experimental::DetectorVolume> constructGapVolume(
    const std::shared_ptr<Acts::Experimental::DetectorVolume>& vol1,
    const std::shared_ptr<Acts::Experimental::DetectorVolume>& vol2,
    std::size_t primaryIdx, std::size_t longIdx, std::size_t shortIdx,
    std::size_t id, const Acts::GeometryContext& gctx) {
  Acts::Vector3 translation1 = vol1->transform(gctx).translation();
  std::vector<double> bounds1 = vol1->volumeBounds().values();

  Acts::Vector3 translation2 = vol2->transform(gctx).translation();
  std::vector<double> bounds2 = vol2->volumeBounds().values();

  if (translation1[longIdx] != translation2[longIdx] ||
      translation1[shortIdx] != translation2[shortIdx]) {
    throw std::runtime_error(
        "Gap volumes can only be constructed between volumes aligned in "
        "non-primary directions");
  }
  if (bounds1[longIdx] != bounds2[longIdx] ||
      bounds1[shortIdx] != bounds2[shortIdx]) {
    throw std::runtime_error(
        "Gap volumes can only be constructed between volumes of equal bounds "
        "along "
        "non-primary directions");
  }
  if (translation1[primaryIdx] + bounds1[primaryIdx] >=
      translation2[primaryIdx] - bounds2[primaryIdx]) {
    throw std::runtime_error(
        "Volumes " + vol1->name() + " and " + vol2->name() +
        " are touching or overlapping, can't construct gap");
  }

  double primaryCenter = (translation2[primaryIdx] - bounds2.at(primaryIdx) +
                          translation1[primaryIdx] + bounds1.at(primaryIdx)) /
                         2.0;
  double primaryHalf = (translation2[primaryIdx] - bounds2.at(primaryIdx) -
                        translation1[primaryIdx] - bounds1.at(primaryIdx)) /
                       2.0;

  double longCenter = translation1[longIdx];
  double shortCenter = translation1[shortIdx];

  double longHalf = bounds1[longIdx];
  double shortHalf = bounds1[shortIdx];

  return constructVolume(primaryHalf, longHalf, shortHalf, primaryCenter,
                         longCenter, shortCenter, primaryIdx, longIdx, shortIdx,
                         "gapVol", id, {}, gctx);
};

void constructGaps(
    const Acts::GeometryContext& gctx, std::size_t primaryIdx,
    std::size_t longIdx, std::size_t shortIdx,
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>&
        detectorVolumes,
    std::size_t idPrefix) {
  std::sort(detectorVolumes.begin(), detectorVolumes.end(),
            [&](const auto& vol1, const auto& vol2) {
              return (vol1->transform(gctx).translation()[primaryIdx] <
                      vol2->transform(gctx).translation()[primaryIdx]);
            });
  std::size_t gapCounter = 0;
  std::size_t nGaps = detectorVolumes.size() - 1;
  for (std::size_t i = 0; i < nGaps; i++) {
    const auto& vol1 = detectorVolumes.at(i);
    const auto& vol2 = detectorVolumes.at(i + 1);

    std::size_t id = idPrefix + gapCounter;
    detectorVolumes.push_back(constructGapVolume(vol1, vol2, primaryIdx,
                                                 longIdx, shortIdx, id, gctx));
    gapCounter++;
  }
  std::sort(detectorVolumes.begin(), detectorVolumes.end(),
            [&](const auto& vol1, const auto& vol2) {
              return (vol1->transform(gctx).translation()[primaryIdx] <
                      vol2->transform(gctx).translation()[primaryIdx]);
            });
}

std::shared_ptr<Acts::PlaneSurface> constructSurface(
    const SurfaceParameters& pars,
    const std::shared_ptr<Acts::RectangleBounds>& bounds,
    const Acts::Vector3& primaryDir, const Acts::Vector3& longDir,
    const Acts::Vector3& shortDir) {
  Acts::Transform3 surfTransform = Acts::Transform3::Identity();

  Acts::RotationMatrix3 surfToWorldRotationPrimary =
      Acts::AngleAxis3(pars.toWorldAnglePrimary, primaryDir).toRotationMatrix();
  Acts::RotationMatrix3 surfToWorldRotationLong =
      Acts::AngleAxis3(pars.toWorldAngleLong, longDir).toRotationMatrix();
  Acts::RotationMatrix3 surfToWorldRotationShort =
      Acts::AngleAxis3(pars.toWorldAngleShort, shortDir).toRotationMatrix();

  surfTransform.translate(pars.toWorldTranslation);

  surfTransform.rotate(surfToWorldRotationPrimary);
  surfTransform.rotate(surfToWorldRotationLong);
  surfTransform.rotate(surfToWorldRotationShort);

  return Acts::Surface::makeShared<Acts::PlaneSurface>(surfTransform, bounds);
};
