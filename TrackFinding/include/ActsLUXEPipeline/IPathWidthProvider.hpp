#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"

/// @brief Interface class for providing path width
/// for a given surface to perform the path seeding
class IPathWidthProvider {
    public:
        virtual ~IPathWidthProvider() = default;

        virtual std::pair<Acts::ActsScalar,Acts::ActsScalar>
        getPathWidth(
            const Acts::GeometryContext& gctx, 
            Acts::GeometryIdentifier geoId) const = 0;
};
