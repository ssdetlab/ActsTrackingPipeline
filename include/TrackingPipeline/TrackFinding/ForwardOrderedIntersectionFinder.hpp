#pragma once

#include "Acts/Surfaces/Surface.hpp" 
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

class ForwardOrderedIntersectionFinder {
    public:
        Acts::ActsScalar m_tol = 1e-4;

        std::vector<const Acts::Surface*> m_surfaces;

        std::vector<std::pair<Acts::GeometryIdentifier, Acts::Vector3>> operator()(
            const Acts::GeometryContext& geoCtx, const Acts::Vector3& position,
            const Acts::Vector3& direction, [[maybe_unused]] const Acts::ActsScalar& Pmag = 0,
            [[maybe_unused]] const Acts::ActsScalar& charge = 0) const {
                std::vector<std::pair<Acts::GeometryIdentifier, Acts::Vector3>> sIntersections;

                // Intersect the surfaces
                for (auto& surface : m_surfaces) {
                    // Get the intersection
                    auto sMultiIntersection = surface->intersect(
                        geoCtx, position, direction,
                        Acts::BoundaryTolerance::Infinite());

                    // Take the closest
                    auto closestForward = sMultiIntersection.closestForward();

                    // Store if the intersection is reachable
                    if (closestForward.status() == Acts::IntersectionStatus::reachable &&
                        closestForward.pathLength() > 0.0) {
                            sIntersections.push_back(
                                {closestForward.object()->geometryId(), closestForward.position()});
                            continue;
                    }
                }

                return sIntersections;
        }
};
