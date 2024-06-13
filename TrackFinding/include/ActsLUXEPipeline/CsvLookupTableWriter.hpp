#pragma once

#include "ActsLUXEPipeline/IWriter.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp" 
#include "ActsLUXEPipeline/DataContainers.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/EventData/SourceLink.hpp"

#include <vector>
#include <fstream>

class CsvLookupTableWriter : public IWriter {
    using eAxis = Acts::detail::EquidistantAxis;
    using eGrid = Acts::Grid<std::vector<Acts::ActsScalar>, eAxis>;

    public:
        struct BinningParameters {
            int nBins;
            double min;
            double max;
        };

        /// @brief The nested configuration struct
        struct Config {
            /// The input collection
            std::string inputCollection = "Measurements";
            /// Surface accessor
            Acts::SourceLinkSurfaceAccessor surfaceAccessor;
            /// First layer extent
            Acts::Extent firstLayerExtent;
            /// Last layer extent
            Acts::Extent lastLayerExtent;
            /// The output file name
            std::string filePath = "";
            /// The binning parameters
            BinningParameters XFirst;
            BinningParameters ZFirst;
        };

        CsvLookupTableWriter() = delete;

        /// @brief Constructor
        CsvLookupTableWriter(const Config& config, Acts::Logging::Level level)
            : m_cfg(config) {
                // Construct a binned grid for each 
                // parameter pair
                eAxis XFirstAxisE(
                    m_cfg.XFirst.min, 
                    m_cfg.XFirst.max, 
                    m_cfg.XFirst.nBins);
                m_EXFirstGrid = std::make_shared<eGrid>(XFirstAxisE);

                eAxis XFirstAxisX(
                    m_cfg.XFirst.min, 
                    m_cfg.XFirst.max, 
                    m_cfg.XFirst.nBins);
                m_XFirstXLastGrid = std::make_shared<eGrid>(XFirstAxisX);

                eAxis XFirstAxisY(
                    m_cfg.XFirst.min, 
                    m_cfg.XFirst.max, 
                    m_cfg.XFirst.nBins);
                m_XFirstYLastGrid = std::make_shared<eGrid>(XFirstAxisY);

                eAxis ZFirstAxis(
                    m_cfg.ZFirst.min, 
                    m_cfg.ZFirst.max, 
                    m_cfg.ZFirst.nBins);
                m_ZFirstZLastGrid = std::make_shared<eGrid>(ZFirstAxis);

                m_inputMeasurements.initialize(m_cfg.inputCollection);
        }

        ~CsvLookupTableWriter() {
            storeLookUp(m_EXFirstGrid, m_cfg.filePath + "XFirstE_lookup_table.csv");
            storeLookUp(m_XFirstXLastGrid, m_cfg.filePath + "XFirstXLast_lookup_table.csv");
            storeLookUp(m_XFirstYLastGrid, m_cfg.filePath + "XFirstYLast_lookup_table.csv");
            storeLookUp(m_ZFirstZLastGrid, m_cfg.filePath + "ZFirstZLast_lookup_table.csv");
        };

        ProcessCode finalize(const AlgorithmContext &ctx) const {
            return ProcessCode::SUCCESS;
        }

        std::string name() const override { return "CsvLookupTableWriter"; };

        void storeLookUp(std::shared_ptr<eGrid> grid, std::string filename) const {
            std::unordered_map<Acts::ActsScalar, Acts::ActsScalar> lookupTable;
            auto nBins = grid->numLocalBins().at(0);
            std::cout << "Number of bins: " << nBins << std::endl;
            // Loop over bins to build the lookup table
            for (std::size_t binx = 0; binx < nBins; binx++) {
                std::vector<Acts::ActsScalar> values = grid->atLocalBins({binx});
                if (values.empty()) {
                    continue;
                }
                std::sort(values.begin(), values.end());
                Acts::ActsScalar median = values[values.size()/2];
                Acts::ActsScalar x = grid->binCenter({binx}).at(0);
                
                lookupTable[x] = median;
            }
            std::fstream file;
            file.open(filename, std::ios::out);
            for (const auto& entry : lookupTable) {
                file << entry.first << "," << entry.second << "\n";
            }
            file.close();
        }

        /// @brief The execute method        
        ProcessCode write(const AlgorithmContext &ctx) override {
            // Get the input measurements
            // from the context
            auto input = m_inputMeasurements(ctx);

            if (input.empty()) {
                return ProcessCode::SUCCESS;
            }

            std::sort(input.begin(), input.end(), 
                [&](const SimMeasurement& a, const SimMeasurement& b) {
                    Acts::Vector2 aHitLoc(
                        a.truthParameters[Acts::eBoundLoc0],
                        a.truthParameters[Acts::eBoundLoc1]);

                    Acts::Vector2 bHitLoc(
                        b.truthParameters[Acts::eBoundLoc0],
                        b.truthParameters[Acts::eBoundLoc1]);

                    auto aHitGlob = 
                        m_cfg.surfaceAccessor(a.sourceLink)->localToGlobal(
                            ctx.geoContext, 
                            aHitLoc, 
                            Acts::Vector3(0, 1, 0));
                    
                    auto bHitGlob =
                        m_cfg.surfaceAccessor(b.sourceLink)->localToGlobal(
                            ctx.geoContext, 
                            bHitLoc, 
                            Acts::Vector3(0, 1, 0));

                    return aHitGlob.y() < bHitGlob.y();
                });

            auto firstHit = input.front();
            Acts::Vector2 locFirstHit(
                firstHit.truthParameters[Acts::eBoundLoc0],
                firstHit.truthParameters[Acts::eBoundLoc1]);

            auto globFirstHit = 
                m_cfg.surfaceAccessor(firstHit.sourceLink)->localToGlobal(
                    ctx.geoContext, 
                    locFirstHit, 
                    Acts::Vector3(0, 1, 0));

            auto me = 0.511 * Acts::UnitConstants::MeV;
            if (m_cfg.firstLayerExtent.contains(globFirstHit)) {
                auto lastHit = input.back();

                Acts::Vector2 locLastHit(
                    lastHit.truthParameters[Acts::eBoundLoc0],
                    lastHit.truthParameters[Acts::eBoundLoc1]);

                auto globLastHit = 
                    m_cfg.surfaceAccessor(lastHit.sourceLink)->localToGlobal(
                        ctx.geoContext, 
                        locLastHit, 
                        Acts::Vector3(0, 1, 0));

                Acts::ActsScalar E = std::hypot(
                    1/firstHit.truthParameters[Acts::eBoundQOverP], me);

                m_EXFirstGrid->atPosition(
                    std::array<double,1>{
                        globFirstHit.x()}).push_back(E);
                m_XFirstXLastGrid->atPosition(
                    std::array<double,1>{
                        globFirstHit.x()}).push_back(globLastHit.x());
                m_XFirstYLastGrid->atPosition(
                    std::array<double,1>{
                        globFirstHit.x()}).push_back(globLastHit.y());
                m_ZFirstZLastGrid->atPosition(
                    std::array<double,1>{
                        globFirstHit.z()}).push_back(globLastHit.z());
            }

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        Config m_cfg;

        /// The output file
        std::fstream m_file;

        /// Private access to the logging instance
        const Acts::Logger &logger() const { return *m_logger; }

        ReadDataHandle<SimMeasurements> m_inputMeasurements
            {this, "InputMeasurements"};

        /// Histograms to handle the binning
        std::shared_ptr<eGrid> m_EXFirstGrid;
        std::shared_ptr<eGrid> m_XFirstXLastGrid;
        std::shared_ptr<eGrid> m_XFirstYLastGrid;
        std::shared_ptr<eGrid> m_ZFirstZLastGrid;

        std::unique_ptr<const Acts::Logger> m_logger;
};
