#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"

#include <memory>

#include <unistd.h>

#include "TrackingPipeline/Alignment/AlignmentContext.hpp"

namespace detail {

std::shared_ptr<AlignmentContext::AlignmentStore> makeAlignmentStore(
    const Acts::GeometryContext& gctx,
    const Acts::Experimental::Detector* detector);

std::shared_ptr<AlignmentContext::AlignmentStore> makeAlignmentStore(
    const Acts::GeometryContext& gctx,
    const Acts::Experimental::Detector* detector,
    const Acts::Vector3& globalShift,
    const std::unordered_map<int, Acts::Vector3>& localShifts,
    const Acts::Vector3& globalAngles,
    const std::unordered_map<int, Acts::Vector3>& localAngles);

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
    const std::unordered_map<int, Acts::Vector3>& localAnglesStdErr);

}  // namespace detail
