#pragma once

#include "ActsLUXEPipeline/ProcessCode.hpp"
#include "ActsLUXEPipeline/IWriter.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp" 
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"

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
class PhoenixTrackWriter : public IWriter {
    public:
        struct Config {
            /// Name of the fitted track collection
            std::string inputTrackCollection;
            /// Name of the output file
            std::string fileName;
        };

        PhoenixTrackWriter(const PhoenixTrackWriter &) = delete;
        PhoenixTrackWriter(const PhoenixTrackWriter &&) = delete;

        /// Constructor
        ///
        /// @param config The configuration struct of the writer
        /// @param level The log level
        PhoenixTrackWriter(const Config& config, Acts::Logging::Level level)
            : m_cfg(config),
            m_logger(Acts::getDefaultLogger(name(), level)) {
                if (m_cfg.fileName.empty()) {
                    throw std::invalid_argument("Missing filename");
                }

                m_inputTracks.initialize(m_cfg.inputTrackCollection);
        }

        /// Writer name() method
        std::string name() const { return "PhoenixTrackWriter"; }

        /// Write out the tracks
        ProcessCode write(const AlgorithmContext &ctx) override {
            auto inputTracks = m_inputTracks(ctx);

            nlohmann::json jsonTracks = nlohmann::json::object();

            jsonTracks["E320Event"] = nlohmann::json::object();
            
            jsonTracks["E320Event"]["Hits"] = nlohmann::json::object();
            jsonTracks["E320Event"]["Hits"]["TrackerHits"] = nlohmann::json::array();

            jsonTracks["E320Event"]["Tracks"] = nlohmann::json::object();
            jsonTracks["E320Event"]["Tracks"]["TrackerTracks"] = nlohmann::json::array();

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
                            10, 10, 10};    
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
};
