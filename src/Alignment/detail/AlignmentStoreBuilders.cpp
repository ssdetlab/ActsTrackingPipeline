#include "TrackingPipeline/Alignment/detail/AlignmentStoreBuilders.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <chrono>
#include <random>
#include <vector>

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

namespace {

Acts::Vector3 computeDetectorCOM(
    const Acts::GeometryContext& gctx,
    const std::vector<Acts::GeometryIdentifier>& geoIds,
    const Acts::Experimental::Detector* detector) {
  Acts::Vector3 com = Acts::Vector3::Zero();
  for (const auto& geoId : geoIds) {
    const auto* surf = detector->findSurface(geoId);
    com += surf->center(gctx) - com;
  }
  com /= geoIds.size();
  return com;
}

std::unordered_map<Acts::GeometryIdentifier, Acts::Vector3>
computeDetectorLeverArms(const Acts::GeometryContext& gctx,
                         const std::vector<Acts::GeometryIdentifier>& geoIds,
                         const Acts::Experimental::Detector* detector) {
  Acts::Vector3 com = computeDetectorCOM(gctx, geoIds, detector);
  std::unordered_map<Acts::GeometryIdentifier, Acts::Vector3> leverArms;
  for (const auto& geoId : geoIds) {
    const auto* surf = detector->findSurface(geoId);
    leverArms[geoId] = surf->center(gctx) - com;
  }
  return leverArms;
}

}  // namespace

namespace detail {

std::shared_ptr<AlignmentContext::AlignmentStore> makeAlignmentStore(
    const Acts::GeometryContext& gctx,
    const Acts::Experimental::Detector* detector) {
  auto aStore = std::make_shared<AlignmentContext::AlignmentStore>();
  for (const auto& vol : detector->volumes()) {
    for (const auto& surf : vol->surfaces()) {
      aStore->emplace(surf->geometryId(), surf->transform(gctx));
    }
  }
  return aStore;
}

std::shared_ptr<AlignmentContext::AlignmentStore> makeAlignmentStore(
    const Acts::GeometryContext& gctx,
    const Acts::Experimental::Detector* detector,
    const Acts::Vector3& globalShift,
    const std::unordered_map<int, Acts::Vector3>& localShifts,
    const Acts::Vector3& globalAngles,
    const std::unordered_map<int, Acts::Vector3>& localAngles) {
  std::vector<Acts::GeometryIdentifier> storeIds;
  for (const auto& vol : detector->volumes()) {
    for (const auto& surf : vol->surfaces()) {
      Acts::GeometryIdentifier geoId = surf->geometryId();
      int id = geoId.sensitive();
      if (localShifts.contains(id) && localAngles.contains(id)) {
        storeIds.push_back(geoId);
      }
    }
  }

  Acts::Vector3 com = computeDetectorCOM(gctx, storeIds, detector);
  std::unordered_map<Acts::GeometryIdentifier, Acts::Vector3> leverArms =
      computeDetectorLeverArms(gctx, storeIds, detector);

  auto aStore = std::make_shared<AlignmentContext::AlignmentStore>();
  for (const auto& geoId : storeIds) {
    const auto* surf = detector->findSurface(geoId);
    int id = geoId.sensitive();

    Acts::Transform3 nominal = surf->transform(gctx);

    Acts::Transform3 transform = Acts::Transform3::Identity();

    Acts::RotationMatrix3 rigidRotation =
        Acts::AngleAxis3(globalAngles.z(), Acts::Vector3::UnitZ())
            .toRotationMatrix() *
        Acts::AngleAxis3(globalAngles.y(), Acts::Vector3::UnitY())
            .toRotationMatrix() *
        Acts::AngleAxis3(globalAngles.x(), Acts::Vector3::UnitX())
            .toRotationMatrix();
    Acts::RotationMatrix3 nominalRotation = nominal.rotation();
    Acts::RotationMatrix3 localRotation =
        Acts::AngleAxis3(localAngles.at(id).z(), Acts::Vector3::UnitZ())
            .toRotationMatrix() *
        Acts::AngleAxis3(localAngles.at(id).y(), Acts::Vector3::UnitY())
            .toRotationMatrix() *
        Acts::AngleAxis3(localAngles.at(id).x(), Acts::Vector3::UnitX())
            .toRotationMatrix();

    Acts::Vector3 translation =
        (com + globalShift) +
        rigidRotation * (leverArms.at(geoId) + localShifts.at(id));
    transform.translation() = translation;

    Acts::RotationMatrix3 rotation =
        rigidRotation * nominalRotation * localRotation;
    transform.rotate(rotation);

    aStore->emplace(geoId, transform);
  }
  return aStore;
}

std::shared_ptr<AlignmentContext::AlignmentStore> makeAlignmentStore(
    const Acts::GeometryContext& gctx,
    const Acts::Experimental::Detector* detector,
    const Acts::Vector3& globalShiftMean,
    const Acts::Vector3& globalShiftStdErr,
    const std::unordered_map<int, Acts::Vector3>& localShiftsMean,
    const std::unordered_map<int, Acts::Vector3>& localShiftsStdErr,
    const Acts::Vector3& globalAnglesMean,
    const Acts::Vector3& globalAnglesStdErr,
    const std::unordered_map<int, Acts::Vector3>& localAnglesMean,
    const std::unordered_map<int, Acts::Vector3>& localAnglesStdErr) {
  RandomNumbers::Config rnCfg;
  rnCfg.seed =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();

  RandomNumbers randomNumberSvc(rnCfg);
  RandomEngine rng = randomNumberSvc.spawnGenerator();
  std::normal_distribution<> normal(0, 1);

  std::vector<Acts::GeometryIdentifier> storeIds;
  for (const auto& vol : detector->volumes()) {
    for (const auto& surf : vol->surfaces()) {
      Acts::GeometryIdentifier geoId = surf->geometryId();
      int id = geoId.sensitive();
      if (localShiftsMean.contains(id) && localShiftsStdErr.contains(id) &&
          localAnglesMean.contains(id) && localShiftsStdErr.contains(id)) {
        storeIds.push_back(geoId);
      }
    }
  }

  Acts::Vector3 com = computeDetectorCOM(gctx, storeIds, detector);
  std::unordered_map<Acts::GeometryIdentifier, Acts::Vector3> leverArms =
      computeDetectorLeverArms(gctx, storeIds, detector);

  Acts::Vector3 globalShift = globalShiftMean + globalShiftStdErr * normal(rng);

  auto aStore = std::make_shared<AlignmentContext::AlignmentStore>();
  for (const auto& geoId : storeIds) {
    const auto* surf = detector->findSurface(geoId);
    int id = geoId.sensitive();
    Acts::Transform3 nominal = surf->transform(gctx);

    Acts::Transform3 transform = Acts::Transform3::Identity();

    Acts::RotationMatrix3 rigidRotation =
        Acts::AngleAxis3(
            globalAnglesMean.z() + globalAnglesStdErr.z() * normal(rng),
            Acts::Vector3::UnitZ())
            .toRotationMatrix() *
        Acts::AngleAxis3(
            globalAnglesMean.y() + globalAnglesStdErr.y() * normal(rng),
            Acts::Vector3::UnitY())
            .toRotationMatrix() *
        Acts::AngleAxis3(
            globalAnglesMean.x() + globalAnglesStdErr.x() * normal(rng),
            Acts::Vector3::UnitX())
            .toRotationMatrix();
    Acts::RotationMatrix3 nominalRotation = nominal.rotation();
    Acts::RotationMatrix3 localRotation =
        Acts::AngleAxis3(localAnglesMean.at(id).z() +
                             localAnglesStdErr.at(id).z() * normal(rng),
                         Acts::Vector3::UnitZ())
            .toRotationMatrix() *
        Acts::AngleAxis3(localAnglesMean.at(id).y() +
                             localAnglesStdErr.at(id).y() * normal(rng),
                         Acts::Vector3::UnitY())
            .toRotationMatrix() *
        Acts::AngleAxis3(localAnglesMean.at(id).x() +
                             localAnglesStdErr.at(id).x() * normal(rng),
                         Acts::Vector3::UnitX())
            .toRotationMatrix();

    Acts::Vector3 translation =
        (com + globalShift) +
        rigidRotation * (leverArms.at(geoId) + localShiftsMean.at(id) +
                         localShiftsStdErr.at(id) * normal(rng));
    transform.translation() = translation;

    Acts::RotationMatrix3 rotation =
        rigidRotation * nominalRotation * localRotation;
    transform.rotate(rotation);

    aStore->emplace(surf->geometryId(), transform);
  }
  return aStore;
}

}  // namespace detail
