#pragma once 

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"

#include <vector>

namespace LUXEMeasurement {

/// @brief Measurement for the LUXE simulation
/// storing the source link and the truth parameters
struct SimMeasurement {
    /// Source link to be used in the
    /// subsequent algorithms
    Acts::SourceLink sourceLink;
    /// The truth parameters 
    Acts::BoundVector truthParameters;
    /// The true track Ids
    std::int32_t trackId;
};

/// @brief A collection of SimMeasurements
using SimMeasurements = std::vector<SimMeasurement>;

} // namespace LUXEMeasurement
