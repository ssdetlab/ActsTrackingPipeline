#pragma once

#include "ActsLUXEPipeline/SimpleSourceLink.hpp"
#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/LUXEEffectiveMaterial.hpp"
#include "ActsLUXEPipeline/SimMeasurementBoostIO.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"
#include "ActsLUXEPipeline/Generators.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp" 
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"

#include "ActsFatras/Physics/ElectroMagnetic/Scattering.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/BetheBloch.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/BetheHeitler.hpp"

#include <memory>
#include <random>
#include <vector>
#include <iostream>
#include <fstream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

using namespace Acts::UnitLiterals;

struct MeasurementsCreatorAction {
    using result_type = SimMeasurements;
    std::size_t sourceId = 0;

    /// @brief Operator that is callable by an ActionList. The function
    /// collects the surfaces
    ///
    /// @tparam propagator_state_t Type of the propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param [out] result Vector of matching surfaces
    /// @param [in] state State of the propagator
    template<typename propagator_state_t, typename stepper_t,
        typename navigator_t>
    void operator()(
        propagator_state_t &state, const stepper_t &stepper,
        const navigator_t &navigator, result_type &result,
        const Acts::Logger &logger) const {
            // only generate measurements on surfaces
            if (!navigator.currentSurface(state.navigation)) {
                return;
            }
            const Acts::Surface &surface = *navigator.currentSurface(state.navigation);
            const Acts::GeometryIdentifier geoId = surface.geometryId();

            // only generate measurements on sensitive surface
            if (!geoId.sensitive()) {
                ACTS_VERBOSE("Create no measurements on non-sensitive surface " << geoId);
                return;
            }

            // Apply global to local
            Acts::Vector3 pos3 = stepper.position(state.stepping);
            Acts::Vector2 loc =
                surface.globalToLocal(
                    state.geoContext, pos3,
                    stepper.direction(state.stepping)).value();

            // The truth info
            Acts::BoundVector parameters = Acts::BoundVector::Zero();
            parameters[Acts::eBoundLoc0] = loc[Acts::eBoundLoc0];
            parameters[Acts::eBoundLoc1] = loc[Acts::eBoundLoc1];
            const auto &direction = stepper.position(state.stepping);
            parameters[Acts::eBoundPhi] = Acts::VectorHelpers::phi(direction);
            parameters[Acts::eBoundTheta] = Acts::VectorHelpers::theta(direction);
            parameters[Acts::eBoundQOverP] = state.stepping.pars[Acts::eFreeQOverP];
            parameters[Acts::eBoundTime] = state.stepping.pars[Acts::eFreeTime];

            // Construct a particle object
            Acts::Vector4 particlePos = {
                pos3.x(), pos3.y(), pos3.z(), parameters[Acts::eBoundTime]};

            ActsFatras::Particle tempParticle =
                makeParticle(
                    stepper.particleHypothesis(state.stepping).absolutePdg(),
                    stepper.direction(state.stepping),
                    stepper.charge(state.stepping)/parameters[Acts::eBoundQOverP],
                    particlePos);

            // Retrieve the material
            Acts::MaterialSlab material;
            if (!surface.surfaceMaterial()) {
                material = makeSiliconSlab();
            }
            auto bsm = dynamic_cast<const Acts::BinnedSurfaceMaterial*>(surface.surfaceMaterial());
            if (bsm) {
                material = bsm->materialSlab(loc);
            }
            auto hsm = dynamic_cast<const Acts::HomogeneousSurfaceMaterial*>(surface.surfaceMaterial());
            if (hsm) {
                material = hsm->materialSlab(loc);
            }

            // Create the scattering and energy loss processes
            std::random_device rd;
            std::mt19937 gen(rd());
            auto scattering = ActsFatras::GaussianMixtureScattering();
            ActsFatras::BetheBloch BBProcess;
            ActsFatras::BetheHeitler BHProcess;

            // Apply the scattering and energy loss processes
            Acts::BoundVector scatteredParameters = parameters;
            ActsFatras::Particle beforeParticle = tempParticle;
            auto scatter = scattering(gen, material, tempParticle);
            auto BBLoss = BBProcess(gen, material, tempParticle);
            auto BHLoss = BHProcess(gen, material, tempParticle);
            scatteredParameters[Acts::eBoundQOverP] =
                tempParticle.charge()/tempParticle.absoluteMomentum();

            // Create the source link
            Acts::SquareMatrix2 cov = Acts::SquareMatrix2::Identity();
            SimpleSourceLink simpleSL(loc, cov, geoId, sourceId);
            
            // Reset the state
            Acts::BoundSquareMatrix resetCov = Acts::BoundSquareMatrix::Identity();
            scatteredParameters[Acts::eBoundPhi] = tempParticle.phi();
            scatteredParameters[Acts::eBoundTheta] = tempParticle.theta();
            stepper.resetState(state.stepping, scatteredParameters, resetCov, surface);
            
            // Dummy IP parameters
            Acts::CurvilinearTrackParameters ipParameters(
                Acts::Vector4::Zero(),
                Acts::Vector3::Zero(),
                1_e / 1_GeV,
                Acts::BoundSquareMatrix::Identity(),
                Acts::ParticleHypothesis::electron());

            // Create the measurement
            Acts::SourceLink sl(simpleSL);
            SimMeasurement sm{
                sl, 
                scatteredParameters,
                ipParameters,
                static_cast<int>(sourceId)};
            result.push_back(sm);
    }
};

/// Propagate the track create smeared measurements from local coordinates.
class MeasurementsCreator : public IAlgorithm {
    public:
        using Propagator = Acts::Propagator<
            Acts::EigenStepper<>, 
            Acts::Experimental::DetectorNavigator>;
    
        using TrackParameters = Acts::CurvilinearTrackParameters;

        /// @brief The nested configuration struct
        struct Config {
            /// The output collection
            std::string outputCollection = "Measurements";
            /// Vertex generator
            std::shared_ptr<IVertexGenerator> vertexGenerator;
            /// Momentum generator
            std::shared_ptr<IMomentumGenerator> momentumGenerator;
            /// Random number generator
            std::shared_ptr<RandomNumbers> randomNumberSvc;
            /// Number of tracks to be generated
            std::size_t nTracks = 100;
        };

        /// @brief Constructor
        MeasurementsCreator(
            Propagator propagator, 
            const Config& config, 
            Acts::Logging::Level level);
        
        ~MeasurementsCreator() = default;
    
        SimMeasurements createMeasurements(
            const Propagator& propagator,
            const AlgorithmContext& ctx,
            const TrackParameters& trackParameters,
            std::size_t id) const;
    
        /// @brief The execute method
        ProcessCode execute(const AlgorithmContext& ctx) const override;
    
        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;
    
        Propagator m_propagator;
    
        WriteDataHandle<SimMeasurements> m_outputMeasurements
            {this, "OutputMeasurements"};
};
