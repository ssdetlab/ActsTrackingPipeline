#pragma once

#include "TrackingPipeline/TrackFinding/detail/TrackLookupTable.hpp"

#include <map>

/// @brief Class to accumulate and average track lookup tables
///
/// @tparam Grid type for track parameters accumulation
///
/// This class is used to accumulate track parameters in
/// reference layer grids and average them to create a lookup
/// table for track parameter estimation in seeding
///
/// @note Geometry context is left to be handled by the user
/// outside of accumulation
class TrackLookupAccumulator {
    public:
        /// @brief Constructor
        explicit TrackLookupAccumulator(TrackLookupGrid grid);

        /// @brief Add track parameters to the accumulator
        ///
        /// @param ipTrackParameters Track parameters at the IP
        /// @param refTrackParameters Track parameters at the reference layer
        /// @param position Local position of the track hit on the reference layer
        void addTrack(
            const Acts::CurvilinearTrackParameters& ipTrackParameters,
            const Acts::CurvilinearTrackParameters& refTrackParameters,
            const Acts::Vector2& position); 

        /// @brief Finalize the lookup table
        ///
        /// @return Grid with the bin track parameters averaged
        TrackLookupGrid finalizeLookup(); 
    
    private:
        /// @brief Add two track parameters
        ///
        /// @param a First track parameter in the sum
        /// @param b Second track parameter in the sum
        ///
        /// @return Sum of track parameters a + b
        ///
        /// @note Covariances of the track parameters
        /// are not added and instead assumed to be
        /// generated by the same random process for
        /// both a and b, making its averaging redundant
        Acts::CurvilinearTrackParameters addTrackParameters(
            const Acts::CurvilinearTrackParameters& a,
            const Acts::CurvilinearTrackParameters& b); 

        /// Grids to accumulate IP and reference
        /// layer track parameters
        TrackLookupGrid m_grid;
        
        /// Mutex for protecting grid access
        std::mutex m_gridMutex;
        
        /// Map to keep the accumulation count
        /// in the occupied grid bins
        std::map<std::array<std::size_t, TrackLookupGrid::DIM>, std::size_t> m_countGrid;
};
