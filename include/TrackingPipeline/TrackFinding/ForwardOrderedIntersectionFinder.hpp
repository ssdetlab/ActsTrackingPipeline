#pragma once

#include "Acts/Surfaces/Surface.hpp" 
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/EventData/TrackParameters.hpp"

/// @brief Class that finds intersection points 
/// based on the pivot paremeters
///
/// Class takes the provided surfaces and finds 
/// intersection points for each pivot source link
/// based on the reference layer track parameters. 
/// The class is designed to work with layer-representing 
/// surfaces combining the senstive surfaces of the detector
/// on the same binning direction, e.g. z-axis.
class ForwardOrderedIntersectionFinder {
    public:
        /// @brief Nested configuration struct
        struct Config {
            /// Merged layers to base the grid on
            std::vector<const Acts::Surface*> layers;
            /// Tolerance
            Acts::ActsScalar tol = 1e-4;
        };

        /// @brief Construct
        ForwardOrderedIntersectionFinder(const Config& config);

        /// @brief Find track intersections with provided 
        /// surfaces
        std::vector<std::pair<Acts::GeometryIdentifier, Acts::Vector2>> operator()(
            const Acts::GeometryContext&, 
            const Acts::CurvilinearTrackParameters&) const;

    private:
        /// Configuration
        Config m_cfg;
};
