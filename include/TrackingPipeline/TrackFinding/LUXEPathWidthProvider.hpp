#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include "TrackingPipeline/Geometry/LUXEGeometryConstraints.hpp"

namespace LUXETrackFinding {

class LUXEPathWidthProvider {
 public:
  LUXEPathWidthProvider(
      LUXEGeometry::GeometryOptions gOpt,
      std::map<std::int32_t, std::pair<double, double>>
          widths)
      : m_gOpt(gOpt), m_widths(widths) {}

  LUXEGeometry::GeometryOptions m_gOpt;
  std::map<std::int32_t, std::pair<double, double>>
      m_widths;

  std::pair<double, double> getPathWidth(
      const Acts::GeometryContext& gctx,
      const Acts::GeometryIdentifier geoId) const {
    std::int32_t staveId = static_cast<std::int32_t>(geoId.sensitive() - 1);
    return m_widths.at(staveId);
  }
};

}  // namespace LUXETrackFinding
