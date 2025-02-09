#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

#include "TrackingPipeline/Clustering/IClusterFilter.hpp"

/// @brief Filter that filters clusters
/// based on their position in the x-y plane
class HourglassFilter : public IClusterFilter {
 public:
  Acts::SourceLinkSurfaceAccessor surfaceAccessor;

  /// Diagonal hourglass lines
  /// slopes
  double a1 = 5.4;
  double a2 = -6.0;

  /// Tunnel-bounded vertical
  /// shifts
  double b1 = 110.0;
  double b2 = 110.0;

  /// Center-bounded vertical
  /// shift
  double b0 = 190.0;

  // Tunnel in the center
  double tunnel = 1.5;

  // Bounds
  double y0 = 90.3 - 29.94176 / 2;
  double x0 = (y0 - b1) / a1;
  double x1 = (y0 - b2) / a2;

  /// Filter operator
  ///
  /// @param cluster: the cluster to filter
  ///
  /// @return: true if the point is inside the
  /// hourglass shape, false otherwise
  bool operator()(const Acts::GeometryContext& gctx,
                  const SimCluster& cluster) const override {
    auto globalPos = surfaceAccessor(Acts::SourceLink(cluster.sourceLink))
                         ->localToGlobal(gctx, cluster.sourceLink.parameters(),
                                         Acts::Vector3(0, 1, 0));
    double x = globalPos.x();
    double y = -globalPos.z();

    // Cluster has to be under the "cap"
    bool cond0 = (x < 0) ? (y < a1 * x + b0) : (y < a2 * x + b0);

    // Conditions alternate between the two sides of
    // the hourglass shape
    if (x < -tunnel) {
      bool cond1 = y < a1 * x + b1;
      bool cond2 = y > a2 * x + b2;
      bool cond3 = x > x0;

      return cond0 && (cond1 || cond2) && cond3;
    } else if (x > tunnel) {
      bool cond1 = y > a1 * x + b1;
      bool cond2 = y < a2 * x + b2;
      bool cond3 = x < x1;

      return cond0 && (cond1 || cond2) && cond3;
    } else {
      // Tunnel in the center
      return cond0;
    }
  }
};
