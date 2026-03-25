#pragma once

#include "Acts/Detector/Detector.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"

namespace E320 {

std::shared_ptr<const Acts::Experimental::Detector> buildDetector(
    const Acts::GeometryContext& gctx,
    const std::shared_ptr<Acts::IMaterialDecorator>& materialDecorator);

std::shared_ptr<Acts::MagneticFieldProvider> buildMagField(
    const Acts::GeometryContext& gctx);

}  // namespace E320Geometry
