#pragma once

#include "TrackingPipeline/Geometry/LUXEGeometryConstraints.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"

namespace LUXETrackFinding {

using namespace Acts::UnitLiterals;

class LUXESourceLinkGrid {
    public:
        using AxisType = Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Open>;
        using GridType = Acts::Grid<std::vector<Acts::SourceLink>, AxisType, AxisType>;

        /// @brief The nested configuration struct
        struct Config {
            /// Geometry options
            LUXEGeometry::GeometryOptions gOpt;
            /// The number of bins
            std::pair<int,int> bins;
            /// Surface accessor
            Acts::SourceLinkSurfaceAccessor surfaceAccessor;
        };
        
        /// @brief Constructor
        LUXESourceLinkGrid(const Config& config) : m_cfg(config) {};

        /// @brief Destructor
        ~LUXESourceLinkGrid() = default;

        /// @brief Bin the source links
        void initialize(
            const Acts::GeometryContext& gctx, 
            std::vector<Acts::SourceLink> sourceLinks) {
                std::unordered_map<int,GridType> lookupTable;
    
                // Construct a binned grid for each layer
                for (int i = 0; i < m_cfg.gOpt.staveZ.size(); i++) {
                    double xMin = 
                        (i % 2 == 0) ? 
                        m_cfg.gOpt.chipXEven.at(0) - m_cfg.gOpt.chipSizeX/2 - 1_mm :
                        m_cfg.gOpt.chipXOdd.at(0) - m_cfg.gOpt.chipSizeX/2 - 1_mm;
                    double xMax =
                        (i % 2 == 0) ?
                        m_cfg.gOpt.chipXEven.at(8) + m_cfg.gOpt.chipSizeX/2 + 1_mm :
                        m_cfg.gOpt.chipXOdd.at(8) + m_cfg.gOpt.chipSizeX/2 + 1_mm;
        
                    double yMin = - m_cfg.gOpt.chipSizeY/2;
                    double yMax = m_cfg.gOpt.chipSizeY/2;
        
                    AxisType xAxis(xMin, xMax, m_cfg.bins.first);
                    AxisType yAxis(yMin, yMax, m_cfg.bins.second);
    
                    GridType grid(std::make_tuple(xAxis, yAxis));
                    lookupTable.insert({i, grid});
                }
                // Fill the grid with source links
                for (auto& sl : sourceLinks) {
                    auto ssl = sl.get<SimpleSourceLink>();
                    auto id = ssl.geometryId().sensitive();
                    int layer = static_cast<int>(id/10 - 1);

                    Acts::Vector3 globalPos = m_cfg.surfaceAccessor(
                        sl)->localToGlobal(
                            gctx, 
                            ssl.parameters(), 
                            Acts::Vector3{0, 1, 0});

                    auto bin = lookupTable.at(layer).localBinsFromPosition(
                        Acts::Vector2(globalPos.x(), globalPos.z()));
                    lookupTable.at(layer).atLocalBins(bin).push_back(sl);
                }
                m_lookupTables = lookupTable;
        }

        GridType operator()(const Acts::GeometryIdentifier& geoId) const {
            return m_lookupTables.at(static_cast<int>(geoId.sensitive()-1));
        }

    private:
        /// Configuration
        Config m_cfg;

        /// Lookup table collection
        std::unordered_map<int,GridType> m_lookupTables;

};

} // namespace LUXESourceLinkGrid
