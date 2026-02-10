#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <chrono>
#include <cstddef>
#include <memory>
#include <random>

#include <unistd.h>

#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

namespace detail {

using namespace Acts::UnitLiterals;

std::shared_ptr<AlignmentContext::AlignmentStore> makeAlignmentStore(
    const Acts::Experimental::Detector* detector) {
  std::map<int, Acts::Vector3> shifts{
      {10, Acts::Vector3(0_mm, 7_um, 17_um)},
      {12, Acts::Vector3(0_mm, 50_um, -10_um)},
      {14, Acts::Vector3(0_mm, -50_um, 10_um)},
      {16, Acts::Vector3(0_mm, -35_um, 35_um)},
      {18, Acts::Vector3(0_mm, 35_um, -35_um)},
      {20, Acts::Vector3(0_mm, 40_um - 24_um, -30_um - 24_um)},
      {22, Acts::Vector3(0_mm, 40_um + 48_um, -30_um + 124_um)},
      {24, Acts::Vector3(0_mm, 40_um + 67_um, -30_um - 40_um)},
      {26, Acts::Vector3(0_mm, 40_um - 67_um, -30_um + 40_um)},
      {28, Acts::Vector3(0_mm, 40_um - 24_um, -30_um - 100_um)}};

  std::map<int, double> angles{
      {10, 1e-3_rad},  {12, 2e-3_rad}, {14, -1.5e-3_rad}, {16, -0.5e-3_rad},
      {18, -1e-3_rad}, {20, 4e-3_rad}, {22, 6e-3_rad},    {24, -7e-3_rad},
      {26, -3e-3_rad}, {28, 1e-4_rad}};

  auto aStore = std::make_shared<AlignmentContext::AlignmentStore>();
  Acts::GeometryContext gctx;
  for (const auto& vol : detector->volumes()) {
    for (const auto& surf : vol->surfaces()) {
      int geoId = surf->geometryId().sensitive();
      if (!geoId) {
        continue;
      }
      Acts::Transform3 nominal = surf->transform(gctx);
      Acts::Vector3 deltaTranslation = shifts.at(geoId);

      Acts::RotationMatrix3 deltaRotation =
          Acts::AngleAxis3(angles.at(geoId), Acts::Vector3::UnitZ())
              .toRotationMatrix();

      deltaTranslation = nominal.rotation().inverse() * deltaTranslation;
      nominal.translate(deltaTranslation);
      nominal.rotate(deltaRotation);

      aStore->emplace(surf->geometryId(), nominal);
    }
  }
  return aStore;
}

std::shared_ptr<AlignmentContext::AlignmentStore> makeAlignmentStore(
    const Acts::Experimental::Detector* detector, std::size_t longIdx,
    double sigmaTransLong, std::size_t shortIdx, double sigmaTransShort) {
  RandomNumbers::Config rnCfg;
  rnCfg.seed =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();

  RandomNumbers randomNumberSvc(rnCfg);
  RandomEngine rng = randomNumberSvc.spawnGenerator();
  std::normal_distribution<> normal(0, 1);

  auto aStore = std::make_shared<AlignmentContext::AlignmentStore>();
  Acts::GeometryContext gctx;
  for (const auto& vol : detector->volumes()) {
    for (const auto& surf : vol->surfaces()) {
      if (!surf->geometryId().sensitive()) {
        continue;
      }
      Acts::Transform3 nominal = surf->transform(gctx);
      Acts::Vector3 delta = Acts::Vector3::Zero();
      delta[longIdx] = sigmaTransLong * normal(rng);
      delta[shortIdx] = sigmaTransShort * normal(rng);

      delta = nominal.rotation().inverse() * delta;
      nominal.translate(delta);

      aStore->emplace(surf->geometryId(), nominal);
    }
  }
  return aStore;
}

}  // namespace detail
