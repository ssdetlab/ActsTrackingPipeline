#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Definitions/Algebra.hpp"

#include <vector>

class IIntersectionFinder {
    public:
        std::vector<const Acts::Surface*> m_surfaces;

        virtual ~IIntersectionFinder() = default;

        virtual std::vector<Acts::SurfaceIntersection> 
        findIntersections(
            const Acts::GeometryContext& gctx, 
            const Acts::Vector3& position,
            const Acts::Vector3& direction) const = 0;
};
