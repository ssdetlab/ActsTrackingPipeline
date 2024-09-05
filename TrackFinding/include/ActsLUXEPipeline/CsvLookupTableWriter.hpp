#pragma once

#include "ActsLUXEPipeline/IWriter.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp" 
#include "ActsLUXEPipeline/DataContainers.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/EventData/SourceLink.hpp"

#include <vector>
#include <fstream>

class CsvLookupTableWriter : public IWriter {
    using TrackParameters = std::tuple<
        Acts::ActsScalar, 
        Acts::ActsScalar, 
        Acts::Vector3, 
        Acts::Vector3, 
        Acts::Vector3>;

    using Axis = Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Open>;
    using Grid = Acts::Grid<std::vector<TrackParameters>, Axis, Axis>;

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
            /// The output file name
            std::string filePath = "";
            /// The binning parameters
            BinningParameters XFirst;
            BinningParameters YFirst;
        };

        CsvLookupTableWriter() = delete;

        /// @brief Constructor
        CsvLookupTableWriter(const Config& config, Acts::Logging::Level level)
            : m_cfg(config) {
                // Construct a binned grid for each 
                // parameter pair
                Axis XFirstAxis(
                    m_cfg.XFirst.min, 
                    m_cfg.XFirst.max, 
                    m_cfg.XFirst.nBins);

                Axis YFirstAxis(
                    m_cfg.YFirst.min, 
                    m_cfg.YFirst.max, 
                    m_cfg.YFirst.nBins);

                m_TrackGrid = std::make_shared<Grid>(
                    std::make_tuple(XFirstAxis, YFirstAxis));

                m_inputMeasurements.initialize(m_cfg.inputCollection);
        }

        ~CsvLookupTableWriter() {
            auto mean = [](const std::vector<TrackParameters>& values) {
                TrackParameters meanTrack = {
                    0, 
                    0, 
                    Acts::Vector3(0,0,0), 
                    Acts::Vector3(0,0,0), 
                    Acts::Vector3(0,0,0)};
                
                for (const auto& track : values) {
                    std::get<0>(meanTrack) += std::get<0>(track);
                    std::get<1>(meanTrack) += std::get<1>(track);
                    std::get<2>(meanTrack) += std::get<2>(track);
                    std::get<3>(meanTrack) += std::get<3>(track);
                    std::get<4>(meanTrack) += std::get<4>(track);
                }

                std::get<0>(meanTrack) /= values.size();
                std::get<1>(meanTrack) /= values.size();
                std::get<2>(meanTrack) /= values.size();
                std::get<3>(meanTrack) /= values.size();
                std::get<4>(meanTrack) /= values.size();

                return meanTrack;
            };

            std::map<
                std::pair<Acts::ActsScalar, Acts::ActsScalar>, 
                TrackParameters> lookupTable{{}};
            auto nBinsX = m_TrackGrid->numLocalBins().at(0);
            auto nBinsY = m_TrackGrid->numLocalBins().at(1);
            // Loop over bins to build the lookup table
            for (std::size_t binx = 0; binx < nBinsX; binx++) {
                for (std::size_t biny = 0; biny < nBinsY; biny++) {
                    auto values = m_TrackGrid->atLocalBins({binx,biny});
                    if (values.empty()) {
                        continue;
                    }

                    TrackParameters meanTrack = mean(values);

                    auto [x,y] = m_TrackGrid->binCenter({binx,biny});

                    lookupTable[{x,y}] = meanTrack;
                }
            }
            std::fstream file;
            file.open(m_cfg.filePath, std::ios::out);
            for (const auto& entry : lookupTable) {
                auto [x,y] = entry.first;
                auto [charge, P, vertex, dirIp, dirFTL] = entry.second;
                file << x << "," << y << "," << charge << "," << P << "," 
                    << vertex.x() << "," << vertex.y() << "," << vertex.z() << ","
                    << dirIp.x() << "," << dirIp.y() << "," << dirIp.z() << ","
                    << dirFTL.x() << "," << dirFTL.y() << "," << dirFTL.z() << std::endl;
            }
            file.close();
        };

        ProcessCode finalize(const AlgorithmContext &ctx) const {
            return ProcessCode::SUCCESS;
        }

        std::string name() const override { return "CsvLookupTableWriter"; };

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

                Acts::Vector3 direction = (globLastHit - globFirstHit).normalized();

                m_TrackGrid->atPosition(
                    std::array<double,2>{
                        globFirstHit.x(), globFirstHit.z()}).push_back(
                            {firstHit.ipParameters.charge(), 
                             firstHit.ipParameters.absoluteMomentum(),
                             firstHit.ipParameters.position(),
                             firstHit.ipParameters.direction(), 
                             direction});
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
        std::shared_ptr<Grid> m_TrackGrid;

        std::unique_ptr<const Acts::Logger> m_logger;
};
