#pragma once

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Simulation/IMeasurementGenerator.hpp"
#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/IDigitizer.hpp"

#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"

#include <memory>
#include <vector>

/// @brief Class creating measurements along 
/// simulated tracks
///
/// Algorithm propagates tracks thorugh the geometry and 
/// create measurements on the encountered senstive 
/// surfaces
class MeasurementsCreator : public IMeasurementGenerator {
    public:
        using Propagator = Acts::Propagator<
            Acts::EigenStepper<>, 
            Acts::Experimental::DetectorNavigator>;

        using TrackParameters = Acts::CurvilinearTrackParameters;

        /// @brief Nested configuration struct
        struct Config {
            /// Vertex generator
            std::shared_ptr<IVertexGenerator> vertexGenerator;
            /// Momentum generator
            std::shared_ptr<IMomentumGenerator> momentumGenerator;
            /// Digitizer
            std::shared_ptr<IDigitizer> hitDigitizer;
            /// Maximum number of steps
            /// to propagate 
            std::size_t maxSteps;
            /// Is background flag
            bool isSignal;
        };

        /// @brief Constructor
        MeasurementsCreator(
            const Propagator propagator,
            const Config& config);

        /// @brief Propagate track and create measurements    
        std::tuple<std::vector<Acts::SourceLink>, SimClusters>
        gen(
            const AlgorithmContext& ctx,
            RandomEngine& rng,
            std::size_t id) const override;
    
        /// @brief Readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        /// Configuration
        Config m_cfg;
    
        /// Propagator instance
        Propagator m_propagator;
};
