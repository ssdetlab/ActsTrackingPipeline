#pragma once

#include "ActsLUXEPipeline/E320GeometryConstraints.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"

#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"

namespace E320TrackFinding {

using namespace Acts::UnitLiterals;

class E320SourceLinkGridConstructor {
    public:
        using AxisType = Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Open>;
        using GridType = Acts::Grid<std::vector<Acts::SourceLink>, AxisType, AxisType>;

        /// @brief The nested configuration struct
        struct Config {
            /// Geometry options
            E320Geometry::GeometryOptions gOpt;
            /// The number of bins
            std::pair<int,int> bins;
            /// Surface accessor
            Acts::SourceLinkSurfaceAccessor surfaceAccessor;
        };
        
        /// @brief Constructor
        E320SourceLinkGridConstructor(const Config& config) : m_cfg(config) {};

        /// @brief Destructor
        ~E320SourceLinkGridConstructor() = default;

        /// @brief Bin the source links
        std::unordered_map<Acts::GeometryIdentifier,GridType> constructGrid(
            const Acts::GeometryContext& gctx, 
            std::vector<Acts::SourceLink> sourceLinks) {
                std::unordered_map<Acts::GeometryIdentifier,GridType> lookupTable;
    
                // Construct a binned grid for each layer
                for (int i = 0; i < m_cfg.gOpt.staveZ.size(); i++) {
                    double yMin = 
                        -m_cfg.gOpt.chipY.at(8) - m_cfg.gOpt.chipSizeY/2 - 1_mm;
                    double yMax =
                        -m_cfg.gOpt.chipY.at(0) + m_cfg.gOpt.chipSizeY/2 + 1_mm;
        
                    double xMin = -m_cfg.gOpt.chipSizeX/2;
                    double xMax = m_cfg.gOpt.chipSizeX/2;
        
                    AxisType xAxis(xMin, xMax, m_cfg.bins.first);
                    AxisType yAxis(yMin, yMax, m_cfg.bins.second);
    
                    Acts::GeometryIdentifier layerId;
                    layerId.setSensitive(i + 1);

                    GridType grid(std::make_tuple(xAxis, yAxis));
                    lookupTable.insert({layerId, grid});
                }
                // Fill the grid with source links
                for (auto& sl : sourceLinks) {
                    auto ssl = sl.get<SimpleSourceLink>();
                    auto sourceLinkId = ssl.geometryId().sensitive();

                    Acts::GeometryIdentifier layerId;
                    int layer = static_cast<int>(sourceLinkId/10 - 1);
                    layerId.setSensitive(layer + 1);

                    Acts::Vector3 globalPos = m_cfg.surfaceAccessor(
                        sl)->localToGlobal(
                            gctx, 
                            ssl.parameters, 
                            Acts::Vector3{0, 1, 0});

                    auto bin = lookupTable.at(layerId).localBinsFromPosition(
                        Acts::Vector2(globalPos.x(), globalPos.z()));
                    lookupTable.at(layerId).atLocalBins(bin).push_back(sl);
                }

                return lookupTable;
        }

    private:
        /// Configuration
        Config m_cfg;
};

} // namespace E320SourceLinkGrid
