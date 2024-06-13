#pragma once

#include "ActsLUXEPipeline/IIntersectionFinder.hpp"

class ForwardOrderedIntersectionFinder : public IIntersectionFinder {
    public:
        Acts::ActsScalar tol = 1e-4;

        std::vector<Acts::SurfaceIntersection> 
        findIntersections(
            const Acts::GeometryContext& gctx, 
            const Acts::Vector3& position,
            const Acts::Vector3& direction) const override {
                std::vector<Acts::SurfaceIntersection> sIntersections;
                // Intersect the surfaces
                for (auto& surface : m_surfaces) {

                    // Get the intersection
                    auto sMultiIntersection = 
                        surface->intersect(
                            gctx, 
                            position, 
                            direction,
                            Acts::BoundaryCheck(true), 
                            tol);
                
                    // Take the closest
                    auto closestForward = sMultiIntersection.closestForward();

                    if (closestForward.status() == 
                            Acts::IntersectionStatus::reachable &&
                        closestForward.pathLength() > 0.0) {
                            sIntersections.push_back(closestForward);
                            continue;
                    }
                }
                // Sort the intersection along the pathlength
                std::sort(sIntersections.begin(), sIntersections.end(),
                    &Acts::SurfaceIntersection::pathLengthOrder);
                return sIntersections;
        }
};
