#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/EventDataView3D.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"

class EventVisualizer : public IAlgorithm {
        public:
        /// @brief The nested configuration struct
        struct Config {
            /// The input collections
            std::string inputCollectionSeeds = "Seeds";
            std::string inputCollectionTracks = "Tracks";
            /// Output path
            std::string outputPath = "";
            /// The number of tracks to be visualized
            std::uint32_t nTracks = 10;
            /// Detector
            const Acts::Experimental::Detector* detector;
            /// Visualization options
            bool visualizeVolumes = true;
            bool visualizeHits = true;
            bool visualizeTracks = true;
        };

        /// @brief Constructor
        EventVisualizer(Config config, Acts::Logging::Level level)
            : IAlgorithm("EventVisualizer", level),
            m_cfg(std::move(config)) {
                m_inputSeeds.initialize(m_cfg.inputCollectionSeeds);
                m_inputTracks.initialize(m_cfg.inputCollectionTracks);
        }
        ~EventVisualizer() = default;

        /// @brief The execute method        
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            // Get the input seeds
            // from the context
            auto inputSeeds = m_inputSeeds(ctx);
            auto inputTracks = m_inputTracks(ctx);

            Acts::ObjVisualization3D obj;
            Acts::Vector3 origin(0, 0, 0);
            Acts::GeometryView3D::drawArrowForward(
                obj, origin, Acts::Vector3(1000, 0, 0),
                1000, 10, Acts::ViewConfig({255, 0, 0}));
            Acts::GeometryView3D::drawArrowForward(
                obj, origin, Acts::Vector3(0, 4000, 0),
                1000, 10, Acts::ViewConfig({0, 255, 0}));
            Acts::GeometryView3D::drawArrowForward(
                obj, origin, Acts::Vector3(0, 0, 100),
                1000, 10, Acts::ViewConfig({0, 0, 255}));

            if (m_cfg.visualizeVolumes) {
                for (auto& vol : m_cfg.detector->rootVolumes()) {
                    Acts::GeometryView3D::drawDetectorVolume(obj, *vol, ctx.geoContext);
                }
            }
            else {
                for (auto& vol : m_cfg.detector->rootVolumes()) {
                    for (auto surf : vol->surfaces()) {
                        Acts::GeometryView3D::drawSurface(obj, *surf, ctx.geoContext);
                    }
                }
            }
            int k = 1;
            if (m_cfg.visualizeHits) {
                for (const auto& seed : inputSeeds) {
                    std::cout << "-----------------------------------" << std::endl;
                    std::cout << "Drawing seed of size " << seed.sourceLinks.size() << std::endl;
                    for (const auto& sl : seed.sourceLinks) {
                        auto ssl = sl.get<SimpleSourceLink>();
                        Acts::Vector2 hit = ssl.parameters;
                        Acts::SquareMatrix2 cov = ssl.covariance;
                        cov *= 1e6;
                        auto surf = 
                            *m_cfg.detector->sensitiveHierarchyMap().find(ssl.geometryId());

                        std::cout << "Surface " << surf->center(ctx.geoContext).transpose() << std::endl;

                        Acts::EventDataView3D::drawMeasurement(
                            obj, hit, cov, surf->transform(ctx.geoContext), 
                            1, Acts::ViewConfig({0, 0, 255}));
                    }
                    std::cout << "-----------------------------------" << std::endl;
                    k++;
                    if (k > m_cfg.nTracks) {
                        break;
                    }
                }
            }
            k = 1;
            if (m_cfg.visualizeTracks) {
                std::cout << "Input tracks size " << inputTracks.size() << std::endl;
                for (const auto& track : inputTracks) {
                    std::vector<Acts::Vector3> trackPoints;
                    for (auto state : track.trackStatesReversed()) {
                        Acts::Vector2 hit(state.parameters()[Acts::eBoundLoc0], 
                            state.parameters()[Acts::eBoundLoc1]);
                        Acts::SquareMatrix2 cov = state.covariance().block<2, 2>(0, 0);
                        cov *= 1e6;
                        auto surf = 
                            *m_cfg.detector->sensitiveHierarchyMap().find(state.referenceSurface().geometryId());

                        Acts::EventDataView3D::drawMeasurement(
                            obj, hit, cov, surf->transform(ctx.geoContext), 
                            1, Acts::ViewConfig({0, 255, 0}));                         

 
                        Acts::Vector3 start = 
                            state.referenceSurface().center(ctx.geoContext) +
                            Acts::Vector3(hit[Acts::eBoundLoc0], hit[Acts::eBoundLoc1], 0);

                        trackPoints.push_back(start);
                    }
                    trackPoints.push_back(origin);
                    for (int i = 0; i < trackPoints.size() - 1; i++) {
                        std::cout << "Drawing track segment (" 
                        << trackPoints[i].transpose() 
                        << ") -> ("
                        << trackPoints[i].transpose() 
                        << ")"
                        << std::endl;
                        Acts::GeometryView3D::drawSegment(
                            obj, trackPoints[i], trackPoints[i + 1],
                            Acts::ViewConfig({0, 255, 0}));
                    }

                    k++;
                    if (k > m_cfg.nTracks) {
                        break;
                    }
                }
            }

            obj.write(m_cfg.outputPath);

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;

        ReadDataHandle<Seeds> m_inputSeeds
            {this, "InputSeeds"};

        ReadDataHandle<Acts::TrackContainer<
            Acts::VectorTrackContainer,
            Acts::VectorMultiTrajectory,
            std::shared_ptr>> m_inputTracks
                {this, "InputTracks"};  

};
