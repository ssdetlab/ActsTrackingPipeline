#pragma once

#include "TrackingPipeline/Clustering/IClusterFilter.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/EventData/SourceLink.hpp"

/// @brief Filter that filters clusters 
/// based on their position in the x-y plane
class HourglassFilter : public IClusterFilter {
    public:
        Acts::SourceLinkSurfaceAccessor surfaceAccessor;

        // Horizontal line
        double a0 = 0.0;
        double b0 = 220.0;
    
        // Diagonal hourglass lines
        double a1 = 5.4;
        double b1 = 110.0;
    
        double a2 = -5.4;
        double b2 = 110.0;
    
        // Tunnel in the center
        double tunnel = 2.0;
    
        /// Filter operator
        ///
        /// @param cluster: the cluster to filter
        /// 
        /// @return: true if the point is inside the 
        /// hourglass shape, false otherwise
        bool operator()(
            const Acts::GeometryContext& gctx, 
            const SimCluster& cluster) const override {
                auto globalPos = surfaceAccessor(
                    Acts::SourceLink(cluster.sourceLink))->localToGlobal(
                        gctx, 
                        cluster.sourceLink.parameters(),
                        Acts::Vector3(0, 1, 0)
                    );
                Acts::ActsScalar x = globalPos.x();
                Acts::ActsScalar y = -globalPos.z();
    
                // Filter out the top part of the tracking plane
                bool cond0 = y < a0 * x + b0;
        
                // Conditions alternate between the two sides of
                // the hourglass shape
                if (x < -tunnel) {
                    bool cond1 = y < a1 * x + b1;
                    bool cond2 = y > a2 * x + b2;
        
                    return cond0 && (cond1 || cond2);
                }
                if (x > tunnel) {
                    bool cond1 = y > a1 * x + b1;
                    bool cond2 = y < a2 * x + b2;
        
                    return cond0 && (cond1 || cond2);
                }
                else {
                    // Tunnel in the center
                    return cond0;
                }
        }
};
