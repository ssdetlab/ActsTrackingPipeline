#pragma once

#include "ActsLUXEPipeline/LUXEGeometryConstraints.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

namespace LUXETrackFinding {

class LUXEPathWidthProvider {
    public:
        LUXEPathWidthProvider(
            LUXEGeometry::GeometryOptions gOpt,
            std::map<std::int32_t, 
                std::pair<
                    Acts::ActsScalar,Acts::ActsScalar>> widths)
            : m_gOpt(gOpt), m_widths(widths) {}

        LUXEGeometry::GeometryOptions m_gOpt;
        std::map<std::int32_t, 
            std::pair<
                Acts::ActsScalar,Acts::ActsScalar>> m_widths;

        std::pair<Acts::ActsScalar,Acts::ActsScalar>
        getPathWidth(
            const Acts::GeometryContext& gctx, 
            const Acts::GeometryIdentifier geoId) const {
                std::int32_t staveId = static_cast<std::int32_t>(geoId.sensitive() - 1);
                return m_widths.at(staveId);
        }
};

} // namespace LUXETrackFinding
