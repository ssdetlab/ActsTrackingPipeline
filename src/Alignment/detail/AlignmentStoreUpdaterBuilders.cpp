#include "TrackingPipeline/Alignment/detail/AlignmentStoreUpdaterBuilders.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <unistd.h>

namespace detail {

ActsAlignment::AlignedTransformUpdater makeGlobalAlignmentUpdater(
    AlignmentContext& alignCtx) {
  return [&alignCtx](Acts::DetectorElementBase* element,
                     const Acts::GeometryContext& gctx,
                     const Acts::Vector3& deltaTranslation,
                     const Acts::Vector3& deltaAngles) {
    const auto& surf = element->surface();
    Acts::Transform3 newTransform = Acts::Transform3::Identity();
    const Acts::Transform3& oldTransform = surf.transform(gctx);

    const Acts::Vector3& oldCenter = oldTransform.translation();
    const Acts::RotationMatrix3& oldRotation = oldTransform.rotation();

    Acts::Vector3 newCenter = oldCenter + deltaTranslation;
    newTransform.translation() = newCenter;
    Acts::RotationMatrix3 newRotation =
        Acts::AngleAxis3(deltaAngles(2), Acts::Vector3::UnitZ())
            .toRotationMatrix() *
        Acts::AngleAxis3(deltaAngles(1), Acts::Vector3::UnitY())
            .toRotationMatrix() *
        Acts::AngleAxis3(deltaAngles(0), Acts::Vector3::UnitX())
            .toRotationMatrix() *
        oldRotation;
    newTransform.rotate(newRotation);

    alignCtx.alignmentStore()[element->surface().geometryId()] = newTransform;

    std::cout << "-----------------------------------\n";
    std::cout << "DELTA ANGLES " << deltaAngles.transpose() << "\n";
    std::cout << "DELTA TRANSLATION " << deltaTranslation.transpose() << "\n";
    std::cout << "SURFACE " << surf.geometryId() << "\n";
    std::cout << "CENTER " << surf.center(gctx).transpose() << "\n";
    std::cout << "NORMAL "
              << surf.normal(gctx, surf.center(gctx), Acts::Vector3::UnitX())
                     .transpose()
              << "\n";
    std::cout << "ROTATION \n" << surf.transform(gctx).rotation() << "\n";
    std::cout << "EXTENT\n"
              << surf.polyhedronRepresentation(gctx, 1000).extent() << "\n";
    return true;
  };
};

ActsAlignment::AlignedTransformUpdater makeLocalAlignmentUpdater(
    AlignmentContext& alignCtx) {
  return [&alignCtx](Acts::DetectorElementBase* element,
                     const Acts::GeometryContext& gctx,
                     const Acts::Vector3& deltaTranslation,
                     const Acts::Vector3& deltaAngles) {
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

    alignCtx.alignmentStore()[element->surface().geometryId()] = newTransform;

    std::cout << "-----------------------------------\n";
    std::cout << "DELTA ANGLES " << deltaAngles.transpose() << "\n";
    std::cout << "DELTA TRANSLATION " << deltaTranslation.transpose() << "\n";
    std::cout << "SURFACE " << surf.geometryId() << "\n";
    std::cout << "CENTER " << surf.center(gctx).transpose() << "\n";
    std::cout << "NORMAL "
              << surf.normal(gctx, surf.center(gctx), Acts::Vector3::UnitX())
                     .transpose()
              << "\n";
    std::cout << "ROTATION \n" << surf.transform(gctx).rotation() << "\n";
    std::cout << "EXTENT\n"
              << surf.polyhedronRepresentation(gctx, 1000).extent() << "\n";
    return true;
  };
};

}  // namespace detail
