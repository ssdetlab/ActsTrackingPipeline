#pragma once

#include "TrackingPipeline/TrackFinding/TrackLookupAccumulator.hpp"
#include "TrackingPipeline/Io/ITrackParamsLookupWriter.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp"

#include <memory>

/// @brief Algorithm to estimate track parameters lookup tables
///
/// This algorithm is used to estimate track parameters lookup tables
/// for track parameter estimation in seeding. The algorithm imposes
/// grids onto the reference tracking layers and accumulates track
/// parameters in the grid bins. The track parameters are then averaged
/// to create a lookup table for track parameter estimation in seeding.
class TrackLookupEstimationAlgorithm : public IAlgorithm {
    public:
        /// @brief Nested configuration struct
        struct Config {
            /// Reference tracking layers
            std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
                refLayers;
            /// Binning of the grid to be emposed
            /// onto the reference layers
            std::pair<std::size_t, std::size_t> bins;
            /// Input SimHit container
            std::string inputClusters = "InputClusters";
            /// Track lookup writers
            std::vector<std::shared_ptr<ITrackParamsLookupWriter>>
                trackLookupGridWriters{};
        };

        /// @brief Constructor
        TrackLookupEstimationAlgorithm(
            const Config& config, Acts::Logging::Level level);
        
        /// @brief Execute method
        ProcessCode execute(const AlgorithmContext& ctx) const override;
        
        /// @brief Finalize method
        ProcessCode finalize() override;
        
        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        /// Configuration
        Config m_cfg;

        /// Input data handles
        ReadDataHandle<SimClusters> m_inputClusters{
            this, "InputSimClusters"};
        
        /// Accumulators for the track parameters
        std::unordered_map<
            Acts::GeometryIdentifier,
            std::unique_ptr<TrackLookupAccumulator>> m_accumulators;
};
