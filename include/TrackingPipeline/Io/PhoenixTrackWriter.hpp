#pragma once

#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp" 
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"

#include "nlohmann/json.hpp"

#include <fstream>

/// @class Json writer abiding the Phoenix format
///
/// @brief Writes out KF tracks in the Phoenix format
class PhoenixWriter : public IWriter {
    public:
        struct Config {
            /// Surface accessor
            Acts::SourceLinkSurfaceAccessor surfaceAccessor;
            /// Fitted track collection
            std::string inputTracks;
            /// Input clusters
            std::string inputClusters;
            /// Name of the output file
            std::string fileName;
        };

        PhoenixWriter(const PhoenixWriter &) = delete;
        PhoenixWriter(const PhoenixWriter &&) = delete;

        /// Constructor
        ///
        /// @param config The configuration struct of the writer
        /// @param level The log level
        PhoenixWriter(const Config& config, Acts::Logging::Level level)
            : m_cfg(config),
            m_logger(Acts::getDefaultLogger(name(), level)) {
                if (m_cfg.fileName.empty()) {
                    throw std::invalid_argument("Missing filename");
                }

                m_inputTracks.initialize(m_cfg.inputTracks);
                m_clusters.initialize(m_cfg.inputClusters);
        }

        /// Writer name() method
        std::string name() const { return "PhoenixWriter"; }

        /// Write out the tracks
        ProcessCode write(const AlgorithmContext &ctx) override {
            auto inputTracks = m_inputTracks(ctx);
            auto inputClusters = m_clusters(ctx);

            nlohmann::json jsonTracks = nlohmann::json::object();

            jsonTracks["E320Event"] = nlohmann::json::object();
            
            jsonTracks["E320Event"]["Hits"] = nlohmann::json::object();
            jsonTracks["E320Event"]["Hits"]["TrackerHits"] = nlohmann::json::array();

            jsonTracks["E320Event"]["Tracks"] = nlohmann::json::object();
            jsonTracks["E320Event"]["Tracks"]["TrackerTracks"] = nlohmann::json::array();

            // Iterate over the clusters
            for (auto cluster : inputClusters) {
                nlohmann::json jsonHit = nlohmann::json::object();

                auto ssl = cluster.sourceLink;
                Acts::Vector3 globalPos = m_cfg.surfaceAccessor(
                    Acts::SourceLink(ssl))->localToGlobal(
                    ctx.geoContext, 
                    ssl.parameters(), 
                    Acts::Vector3{0, 1, 0});

                jsonHit["type"] = "Box";
                jsonHit["pos"] = {
                    globalPos.x(),
                    -globalPos.z(), 
                    globalPos.y(),
                    1, 1, 1};
                jsonTracks["E320Event"]["Hits"]["ClusterHits"].push_back(jsonHit);
            }

            // Iterate over the fitted tracks
            for (int idx = 0; idx < inputTracks.size(); idx++) {
                // Get the track object and the track id
                auto [id,track] = inputTracks.getByIndex(idx);

                nlohmann::json jsonTrack = nlohmann::json::object();
                jsonTrack["color"] = "0x00ff00";
                jsonTrack["pos"] = nlohmann::json::array();

                // Iterate over the track states
                for (auto state : track.trackStatesReversed()) {
                    auto recoPars = state.parameters();
                    Acts::Vector2 recoLocalPos = recoPars.head<2>();
                    Acts::Vector3 recoGlobalPos = state.referenceSurface().localToGlobal(
                        ctx.geoContext, recoLocalPos, Acts::Vector3(1, 0, 0));

                    jsonTrack["pos"].push_back(
                        {recoGlobalPos.x(), 
                        -recoGlobalPos.z(), 
                        recoGlobalPos.y()}
                    );

                    // Skip the states without meaningful information
                    if (state.hasUncalibratedSourceLink()) {
                        auto sl = state.getUncalibratedSourceLink();
                        auto ssl = sl.get<SimpleSourceLink>();

                        Acts::Vector2 hitLocalPos = ssl.parameters();
                        Acts::Vector3 hitGlobalPos = state.referenceSurface().localToGlobal(
                            ctx.geoContext, hitLocalPos, Acts::Vector3(1, 0, 0));

                        nlohmann::json jsonHit = nlohmann::json::object();
                        jsonHit["type"] = "Box";
                        jsonHit["pos"] = {
                            hitGlobalPos.x(), 
                            -hitGlobalPos.z(), 
                            hitGlobalPos.y(),
                            1, 1, 1};    
                        jsonTracks["E320Event"]["Hits"]["TrackerHits"].push_back(jsonHit);
                    }
                }
                jsonTracks["E320Event"]["Tracks"]["TrackerTracks"].push_back(jsonTrack);
            }
            // Evoke the converter
            // And write the file(s)
            std::string fileName = m_cfg.fileName + ".json";
            ACTS_VERBOSE("Writing to file: " << fileName);
            std::ofstream ofj(fileName);
            ofj << std::setw(4) << jsonTracks << std::endl;

            return ProcessCode::SUCCESS;
        }

        /// Readonly access to the config
        const Config& config() const { return m_cfg; }

    private:
        const Acts::Logger& logger() const { return *m_logger; }
        
        /// The logger instance
        std::unique_ptr<const Acts::Logger> m_logger{nullptr};
        
        /// The config of the writer
        Config m_cfg;

        ReadDataHandle<Tracks<
            Acts::VectorTrackContainer,
            Acts::VectorMultiTrajectory>>
                m_inputTracks{this, "InputTracks"};  

        ReadDataHandle<SimClusters> m_clusters{this, "Clusters"};
};
