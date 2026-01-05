#pragma once

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "TrackingPipeline/Geometry/detail/SurfaceParameters.hpp"

std::shared_ptr<Acts::Experimental::DetectorVolume> constructVolume(
    double halfPrimary, double halfLong, double halfShort, double centerPrimary,
    double centerLong, double centerShort, std::size_t primaryIdx,
    std::size_t longIdx, std::size_t shortIdx, const std::string& namePrefix,
    std::size_t id, const std::vector<std::shared_ptr<Acts::Surface>>& surfaces,
    const Acts::GeometryContext& gctx);

std::shared_ptr<Acts::Experimental::DetectorVolume> constructGapVolume(
    const std::shared_ptr<Acts::Experimental::DetectorVolume>& vol1,
    const std::shared_ptr<Acts::Experimental::DetectorVolume>& vol2,
    std::size_t primaryIdx, std::size_t longIdx, std::size_t shortIdx,
    std::size_t id, const Acts::GeometryContext& gctx);

std::shared_ptr<Acts::PlaneSurface> constructSurface(
    const SurfaceParameters& pars,
    const std::shared_ptr<Acts::RectangleBounds>& bounds,
    const Acts::Vector3& primaryDir, const Acts::Vector3& longDir,
    const Acts::Vector3& shortDir);
