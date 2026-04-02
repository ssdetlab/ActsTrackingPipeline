#include "TrackingPipeline/Alignment/LocalAlignmentTransformUpdater.hpp"

#include "Acts/Surfaces/Surface.hpp"

LocalAlignmentTransformUpdater::LocalAlignmentTransformUpdater(
    const Config& cfg, Acts::Logging::Level level)
    : m_cfg(cfg),
      m_logger(
          Acts::getDefaultLogger("LocalAlignmentTransformUpdater", level)) {}

bool LocalAlignmentTransformUpdater::updateAlignmentParameters(
    Acts::DetectorElementBase* element, const Acts::GeometryContext& gctx,
    const Acts::Vector3& deltaTranslation,
    const Acts::Vector3& deltaAngles) const {
  const auto& surf = element->surface();
  Acts::Transform3 newTransform = Acts::Transform3::Identity();
  const Acts::Transform3& oldTransform = surf.transform(gctx);

  const Acts::Vector3& oldCenter = oldTransform.translation();
  const Acts::RotationMatrix3& oldRotation = oldTransform.rotation();

  Acts::Vector3 newCenter = oldCenter + deltaTranslation;
  newTransform.translation() = newCenter;
  Acts::RotationMatrix3 newRotation =
      oldRotation *
      Acts::AngleAxis3(deltaAngles(2), Acts::Vector3::UnitZ())
          .toRotationMatrix() *
      Acts::AngleAxis3(deltaAngles(1), Acts::Vector3::UnitY())
          .toRotationMatrix() *
      Acts::AngleAxis3(deltaAngles(0), Acts::Vector3::UnitX())
          .toRotationMatrix();
  newTransform.rotate(newRotation);

  m_cfg.alignmentStore->store[element->surface().geometryId()] = newTransform;

  ACTS_INFO("Delta angles " << deltaAngles.transpose());
  ACTS_INFO("Delta translation " << deltaTranslation.transpose());
  ACTS_INFO("Surface " << surf.geometryId());
  ACTS_INFO("Center " << surf.center(gctx).transpose());
  ACTS_INFO(
      "Normal " << surf.normal(gctx, surf.center(gctx), Acts::Vector3::UnitX())
                       .transpose());
  ACTS_INFO("Rotation \n" << surf.transform(gctx).rotation());
  ACTS_INFO("Extent\n" << surf.polyhedronRepresentation(gctx, 1000).extent());
  return true;
}
