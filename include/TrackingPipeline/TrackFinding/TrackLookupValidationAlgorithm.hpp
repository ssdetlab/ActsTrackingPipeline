#pragma once

#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp"

#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/SourceLink.hpp"

/// @brief Algorithm setting up comparison of estimated 
/// and truth track parameters
///
/// Algorithm iterates over the input clusters with 
/// truth information and for each hit on the reference 
/// layers estimates track parameters at the IP and 
/// the given reference layer
class TrackLookupValidationAlgorithm : public IAlgorithm {
    public:
        /// @brief Delegate to estimate the IP parameters
        /// and the momentum direction at the reference tracking layer
        ///
        /// @arg Geometry context to use
        /// @arg Pivot source link
        ///
        /// @return Pair of the track parameters at the IP and
        /// the reference tracking layer
        using TrackEstimator = Acts::Delegate<
            std::pair<
                Acts::CurvilinearTrackParameters, 
                Acts::CurvilinearTrackParameters>(
                    const Acts::GeometryContext&, const Acts::SourceLink&)>;
        
        /// @brief Nested configuration struct
        struct Config {
            /// Reference tracking layers
            std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
                refLayers;
            /// Track parameters estimator
            TrackEstimator estimator;
            /// Input Measurement collection
            std::string inputClusters = "InputParticles";
            /// Output IpPars collection
            std::string outputIpPars = "OutputIpPars";
            /// Output RefLayerPars collection
            std::string outputRefLayerPars = "OutputRefLayerPars";
            /// Output IpParsEst collection
            std::string outputIpParsEst = "OutputIpParsEst";
            /// Output RefLayerParsEst collection
            std::string outputRefLayerParsEst = "OutputRefLayerParsEst";
        };

        /// @brief Constructor
        TrackLookupValidationAlgorithm(const Config& config, Acts::Logging::Level level);
        
        /// @brief Execute method
        ProcessCode execute(const AlgorithmContext& ctx) const final;
        
        /// @brief Eeadonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        /// Configuratoin
        Config m_cfg;
        
        /// Input clusters with truth parameters
        ReadDataHandle<SimClusters> m_inputClusters{
            this, "InputSimClusters"};
        
        /// Output track parameters of the hits
        WriteDataHandle<std::vector<
            Acts::CurvilinearTrackParameters>> m_outputIpPars{
                this, "OutputIpPars"};
        WriteDataHandle<std::vector<
            Acts::CurvilinearTrackParameters>> m_outputRefLayerPars{
                this, "OutputRefLayerPars"};
        
        /// Output estimated track parameters
        WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
            m_outputIpParsEst{this, "OutputIpParsEst"};
        WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
            m_outputRefLayerParsEst{this, "OutputRefLayerParsEst"};
};
