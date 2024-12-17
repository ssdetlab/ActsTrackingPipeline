#pragma once

#include "TrackingPipeline/TrackFinding/detail/TrackLookupTable.hpp"

/// @brief Interface for writing track parameter lookup tables
class ITrackParamsLookupWriter {
    public:
        /// Virtual Destructor
        virtual ~ITrackParamsLookupWriter() = default;

        /// Writer method
        ///
        /// @param lookup track lookup to write
        virtual void writeLookup(const TrackLookup& lookup) const = 0;
};
