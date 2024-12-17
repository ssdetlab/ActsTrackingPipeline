#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"

#include <memory>
#include <unordered_map>

using TrackLookupPair =
    std::pair<std::shared_ptr<Acts::CurvilinearTrackParameters>,
              std::shared_ptr<Acts::CurvilinearTrackParameters>>;

/// @brief Track parameters lookup table axis used
/// in the track estimation algorithm
using TrackLookupAxis =
    Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Open>;

/// @brief Track parameters lookup table axis generator
/// used in the track estimation algorithm
using TrackLookupAxisGen = Acts::GridAxisGenerators::EqOpenEqOpen;

/// @brief Lookup grid for track parameters estimation
/// in a given layer
using TrackLookupGrid =
    Acts::Grid<TrackLookupPair, TrackLookupAxis, TrackLookupAxis>;

/// @brief Lookup table for track parameters estimation
/// in the track estimation algorithm
using TrackLookup =
    std::unordered_map<Acts::GeometryIdentifier, TrackLookupGrid>;
