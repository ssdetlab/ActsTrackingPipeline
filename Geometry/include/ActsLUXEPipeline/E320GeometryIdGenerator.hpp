#pragma once

#include "Acts/Detector/interface/IGeometryIdGenerator.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "ActsLUXEPipeline/E320GeometryConstraints.hpp"

#include <any>

namespace E320Geometry {

/// @brief The geometry id generator
/// working as per internal E320 conventions
class E320GeometryIdGenerator final : 
    public Acts::Experimental::IGeometryIdGenerator {
        public:
            /// @brief  Nested config struct
            struct Config {
                /// Container mode
                bool containerMode = false;
                /// Container id (if container mode), will not be incremented
                unsigned int containerId = 0u;
                /// Resetting mode
                bool resetSubCounters = true;
                /// Force override existing ids
                bool overrideExistingIds = false;
                /// Geometry options to extract Ids
                const E320Geometry::GeometryOptions& gOpt;
            };

            /// @brief Nested cache struct
            struct Cache {
                /// Cache count of the volume, for non-container mode
                unsigned int volumeCount = 0u;
                /// Cache count of the layer volume, for container mode
                unsigned int layerCount = 0u;
                /// Cache count of the portal surfaces
                unsigned int portalCount = 0u;
                /// Cache count of passive surfaces
                unsigned int passiveCount = 0u;
                /// Cache count of sensitive surfaces
                unsigned int sensitiveCount = 0u;
            };

            /// @brief Constructor with config
            ///
            /// @param cfg is the geometry configuration object
            /// @param mlogger is the logging instance
            E320GeometryIdGenerator(const Config& cfg,
                std::unique_ptr<const Acts::Logger> mlogger = 
                    Acts::getDefaultLogger("E320GeometryIdGenerator", 
                        Acts::Logging::INFO))
                : m_cfg(cfg), m_logger(std::move(mlogger)) {}
        
            ~E320GeometryIdGenerator() override = default;

            /// @brief Interface method to generate a geometry id cache
            /// @return a geometry id cache wrapped in a std::any object
            IGeometryIdGenerator::GeoIdCache generateCache() const final;

            /// @brief Method for assigning a geometry id to a detector volume
            ///
            /// @param cache is the cache object for e.g. object counting
            /// @param dVolume the detector volume to assign the geometry id to
            void assignGeometryId(
                Acts::Experimental::IGeometryIdGenerator::GeoIdCache& cache,
                Acts::Experimental::DetectorVolume& dVolume) const final;

            /// @brief Method for assigning a geometry id to a portal
            ///
            /// @param cache is the cache object for e.g. object counting
            /// @param portal the portal to assign the geometry id to
            void assignGeometryId(
                Acts::Experimental::IGeometryIdGenerator::GeoIdCache& cache,
                Acts::Experimental::Portal& portal) const final;

            /// @brief Method for assigning a geometry id to a surface
            ///
            /// @param cache is the cache object for e.g. object counting
            /// @param surface the surface to assign the geometry id to
            void assignGeometryId(
                Acts::Experimental::IGeometryIdGenerator::GeoIdCache& cache,
                Acts::Surface& surface) const final;

        private:
            /// @brief Helper method to get the volume id from the cache
            ///
            /// @param cache the provided cache
            /// @param incrementLayer if true, the layer counter is incremented
            ///
            /// @return a valid geometry identifier
            Acts::GeometryIdentifier volumeId(
                Cache& cache, bool incrementLayer = true) const;
        
            /// Configuration object
            Config m_cfg;
        
            /// Private access method to the logger
            const Acts::Logger& logger() const { return *m_logger; }
        
            /// logging instance
            std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace E320Geometry
