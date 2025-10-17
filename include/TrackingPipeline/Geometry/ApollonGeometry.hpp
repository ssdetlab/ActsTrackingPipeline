#pragma once

#include "Acts/Detector/Detector.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <memory>

namespace ApollonGeometry {

std::shared_ptr<const Acts::Experimental::Detector> buildDetector(
    const Acts::GeometryContext& gctx, bool insertReferenceSurface = false);

std::shared_ptr<Acts::MagneticFieldProvider> buildMagField(
    const Acts::GeometryContext& gctx);

}  // namespace ApollonGeometry
