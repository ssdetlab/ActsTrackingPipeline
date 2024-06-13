#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"

class ISourceLinkBinner {
    using eAxis = Acts::detail::EquidistantAxis;
    using eGrid = Acts::Grid<std::vector<Acts::SourceLink>, eAxis, eAxis>;

    public:
        /// @brief Virtual destructor
        virtual ~ISourceLinkBinner() = default;

        /// @brief Interface function to bin source links
        virtual void initialize(
            const Acts::GeometryContext& gctx, 
            std::vector<Acts::SourceLink> sourceLinks) = 0;

        virtual eGrid getLookupTable(const Acts::GeometryIdentifier& geoId) const = 0;

};