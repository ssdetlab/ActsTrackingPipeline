#include "TrackingPipeline/TrackFinding/ForwardOrderedIntersectionFinder.hpp"

#include "Acts/Surfaces/Surface.hpp" 
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include <Acts/Definitions/Algebra.hpp>

ForwardOrderedIntersectionFinder::ForwardOrderedIntersectionFinder(
    const Config& config) : m_cfg(config) {};

std::vector<std::pair<Acts::GeometryIdentifier, Acts::Vector2>> 
ForwardOrderedIntersectionFinder::operator()(
    const Acts::GeometryContext& gctx, 
    const Acts::CurvilinearTrackParameters& refParameters) const {
        std::vector<
            std::pair<
                Acts::GeometryIdentifier, 
                Acts::Vector2>> sIntersections;

        Acts::Vector3 position = refParameters.position(gctx);
        Acts::Vector3 direction = refParameters.direction();

        // Intersect the surfaces
        for (auto& layer : m_cfg.layers) {
            // Get the intersection
            auto sMultiIntersection = 
                layer->intersect(
                    gctx, 
                    position, 
                    direction,
                    Acts::BoundaryTolerance::Infinite());

            // Take the closest
            auto closestForward = 
                sMultiIntersection.closestForward();

            // Store if the intersection is reachable
            if (closestForward.status() == Acts::IntersectionStatus::reachable &&
                closestForward.pathLength() > 0.0) {
                    Acts::Vector2 intersectionLocal =
                        closestForward.object()->globalToLocal(
                            gctx,
                            closestForward.position(), 
                            Acts::Vector3(0, 1, 0)).value();

                    sIntersections.push_back(
                        {closestForward.object()->geometryId(), 
                        intersectionLocal});
                    continue;
            }
        }

        return sIntersections;
}
