#pragma once

#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

namespace E320TrackFinding {

class E320PathWidthProvider {
    public:
        E320PathWidthProvider(
            E320Geometry::GeometryOptions gOpt,
            std::map<std::int32_t, 
                std::pair<
                    Acts::ActsScalar,Acts::ActsScalar>> widths)
            : m_gOpt(gOpt), m_widths(widths) {}

        E320Geometry::GeometryOptions m_gOpt;
        std::map<std::int32_t, 
            std::pair<
                Acts::ActsScalar,Acts::ActsScalar>> m_widths;

        std::pair<Acts::ActsScalar,Acts::ActsScalar>
        operator()(
            const Acts::GeometryContext& /*gctx*/, 
            const Acts::GeometryIdentifier& geoId) const {
                std::int32_t staveId = static_cast<std::int32_t>(geoId.sensitive() - 1);
                return m_widths.at(staveId);
        }
};

} // namespace E320TrackFinding
