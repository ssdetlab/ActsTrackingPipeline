#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Definitions/Algebra.hpp"

class ILookupProvider {
    public:
        /// @brief Virtual destructor
        virtual ~ILookupProvider() = default;

        /// @brief Interface function to get the lookup table
        /// for par_1 and par_2 values
        virtual std::unordered_map<
            Acts::ActsScalar,Acts::ActsScalar> getLookup() const = 0;
};