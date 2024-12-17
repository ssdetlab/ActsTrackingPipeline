#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <map>

class LayerPathWidthProvider {
    public:
        LayerPathWidthProvider(
            std::map<std::int32_t, 
                std::pair<
                    Acts::ActsScalar,Acts::ActsScalar>> widths)
            : m_widths(widths) {}

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
