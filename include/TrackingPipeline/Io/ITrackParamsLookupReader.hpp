#pragma once

#include "TrackingPipeline/TrackFinding/detail/TrackLookupTable.hpp"

/// @brief Interface for reading track parameter lookup tables
class ITrackParamsLookupReader {
    public:
        /// Virtual Destructor
        virtual ~ITrackParamsLookupReader() = default;

        /// Reader method
        ///
        /// @param path the path to the file to read
        virtual TrackLookup readLookup(const std::string& path) const = 0;
};
