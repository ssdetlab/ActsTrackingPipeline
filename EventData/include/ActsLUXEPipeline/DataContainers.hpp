#pragma once 

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackProxy.hpp"

#include <vector>

/// @brief Measurement for the LUXE simulation
/// storing the source link and the truth parameters
struct SimMeasurement {
    /// Source link to be used in the
    /// subsequent algorithms
    Acts::SourceLink sourceLink;
    /// The truth parameters 
    Acts::BoundVector truthParameters;
    /// True vertex 
    Acts::Vector4 trueVertex;
    /// The true track Ids
    std::int32_t trackId;
};

/// @brief A collection of SimMeasurements
using SimMeasurements = std::vector<SimMeasurement>;

struct Seed {
    /// Source links related
    /// to the seed measurements
    std::vector<Acts::SourceLink> sourceLinks;
    /// IP parameters
    Acts::CurvilinearTrackParameters ipParameters;
    /// The track Ids
    std::int32_t trackId;
};

using Seeds = std::vector<Seed>;

template <typename container_t, typename trajectory_t>
struct Tracks {
    using TrackContainer =
        Acts::TrackContainer<container_t, trajectory_t, std::shared_ptr>;
        
    using IdContainer = std::vector<std::int32_t>;

    using ConstTrackProxy = 
        Acts::TrackProxy<container_t, trajectory_t, std::shared_ptr, true>;


    TrackContainer tracks;
    IdContainer trackIds;

    std::pair<std::int32_t, ConstTrackProxy> 
    getByIndex(std::int32_t i) {
        return {trackIds.at(i),tracks.getTrack(i)};
    } 

    std::pair<std::int32_t, ConstTrackProxy> 
    getByTrackId(std::int32_t i) {
        auto it = std::find(trackIds.begin(), trackIds.end(), i);
        if (it != trackIds.end()) {
            return {i,tracks.getTrack(std::distance(trackIds.begin(), it))};
        }
        else {
            return {i,nullptr};
        }
    }

    std::pair<std::int32_t, ConstTrackProxy> 
    begin() {
        return {trackIds.at(0),tracks.getTrack(0)};
    }

    std::pair<std::int32_t, ConstTrackProxy> 
    end() {
        return {trackIds.at(trackIds.size()-1),tracks.getTrack(trackIds.size()-1)};
    }

    std::size_t size() {
        return trackIds.size();
    }
};
