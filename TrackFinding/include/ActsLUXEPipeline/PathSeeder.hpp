#pragma once

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"
#include "ActsLUXEPipeline/ISourceLinkBinner.hpp"
#include "ActsLUXEPipeline/ILookupProvider.hpp"
#include "ActsLUXEPipeline/IPathWidthProvider.hpp"
#include "ActsLUXEPipeline/IIntersectionFinder.hpp" 

#include <unordered_map>

class PathSeeder : public IAlgorithm {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// The input collection
            std::string inputCollection = "SourceLink";
            /// The output collection
            std::string outputCollection = "Seed";
            /// Binned SourceLink provider
            std::shared_ptr<ISourceLinkBinner> sourceLinkBinner;
            /// Lookup providers
            std::shared_ptr<ILookupProvider> EXFirstLookupProvider;
            std::shared_ptr<ILookupProvider> XFirstXLastLookupProvider;
            std::shared_ptr<ILookupProvider> XFirstYLastLookupProvider;
            std::shared_ptr<ILookupProvider> ZFirstZLastLookupProvider;
            /// Surface accessor
            Acts::SourceLinkSurfaceAccessor surfaceAccessor;
            /// Intersection finder
            std::shared_ptr<IIntersectionFinder> intersectionFinder;
            /// Path width provider
            std::shared_ptr<IPathWidthProvider> pathWidthProvider;
            /// First layer extent
            Acts::Extent firstLayerExtent;
            /// The minimum number of hits
            /// for a seed to be created
            std::uint32_t minHits = 4;
            /// The maximum number of hits
            /// for a seed to be created
            std::uint32_t maxHits = 4;
        };

        /// @brief Constructor
        PathSeeder(Config config, Acts::Logging::Level level)
            : IAlgorithm("PathSeeder", level),
            m_cfg(std::move(config)) {
                m_inputMeasurements.initialize(m_cfg.inputCollection);
                m_outputSeeds.initialize(m_cfg.outputCollection);
        }
        ~PathSeeder() = default;

        /// @brief The seeding function
        ///
        /// @param ctx The algorithm context
        /// @param input The input data handle
        /// @param output The output data handle
        ProcessCode execute(const AlgorithmContext& ctx) const override;

    private:
        Acts::ActsScalar
        findClosestValue(
            std::unordered_map<
                Acts::ActsScalar,Acts::ActsScalar>& lookupTable, 
            Acts::ActsScalar x) const;

        Config m_cfg;

        ReadDataHandle<SimMeasurements> m_inputMeasurements
            {this, "InputMeasurements"};

        WriteDataHandle<Seeds> m_outputSeeds
            {this, "OutputSeeds"};
};
