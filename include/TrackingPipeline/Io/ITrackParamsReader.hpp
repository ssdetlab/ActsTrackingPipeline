#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include <vector>

/// @brief Interface for reading track parameters 
class ITrackParamsReader {
    public:
        /// Virtual Destructor
        virtual ~ITrackParamsReader() = default;

        /// Reader method
        ///
        /// @param path the path to the file to read
        virtual std::vector<Acts::CurvilinearTrackParameters> 
            read() = 0;
};
