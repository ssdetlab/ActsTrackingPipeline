#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"

/// @brief Interface class for binning source links
/// into a lookup table -- a grid of source links
class ISourceLinkBinner {
    using eAxis = Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Open>;
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