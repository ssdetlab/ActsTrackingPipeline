#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"
#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/LUXEEffectiveMaterial.hpp"
#include "ActsLUXEPipeline/LUXEMeasurement.hpp"

#include "ActsLUXEPipeline/DataHandle.hpp"

#include "ActsLUXEPipeline/LUXEDataContainers.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
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

/// Propagator action to create measurements.
namespace LUXENavigator {

    using namespace Acts::UnitLiterals;

    struct MeasurementsCreatorAction {

        using result_type = LUXEDataContainer::SimMeasurements;
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
        void operator()(propagator_state_t &state, const stepper_t &stepper,
                        const navigator_t &navigator, result_type &result,
                        const Acts::Logger &logger) const {
            using namespace Acts::UnitLiterals;
            Acts::Vector4 trueVertex;
            if (state.steps == 0) {
                trueVertex = {stepper.position(state.stepping)[0],
                              stepper.position(state.stepping)[1],
                              stepper.position(state.stepping)[2],
                              stepper.time(state.stepping)};
            }
            if (state.steps > 100) {
                std::cout<<"Too many steps, discarding run"<<std::endl;
                stepper.update(state.stepping, Acts::Vector3{0,0,0}, Acts::Vector3{0,1,0},
                               .01, 0);
                return;
            }

            // only generate measurements on surfaces
            if (!navigator.currentSurface(state.navigation)) {
//                result.fullTrack.push_back(stepper.position(state.stepping));
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
                    surface.globalToLocal(state.geoContext, pos3,
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
            Acts::Vector4 particlePos = {pos3[0],pos3[1],pos3[2],parameters[Acts::eBoundTime]};
            ActsFatras::Particle tempParticle =
                    makeParticle(
                            stepper.particleHypothesis(state.stepping).absolutePdg(),stepper.direction(state.stepping),
                            stepper.charge(state.stepping)/parameters[Acts::eBoundQOverP],
                            particlePos);
            std::vector<Acts::MaterialSlab> layerMaterials = {makeSiliconSlab()};
            std::random_device rd;
            std::mt19937 gen(rd());
            auto scattering = ActsFatras::GaussianMixtureScattering();
            ActsFatras::BetheBloch BBProcess;
            ActsFatras::BetheHeitler BHProcess;

            Acts::BoundVector scatteredParameters = parameters;
            ActsFatras::Particle beforeParticle = tempParticle;
            for (auto material:layerMaterials) {
                float p1 = tempParticle.absoluteMomentum();
                auto scatter = scattering(gen, material, tempParticle);
                auto BBLoss = BBProcess(gen, material, tempParticle);
                auto BHLoss = BHProcess(gen, material, tempParticle);
                float p2 = tempParticle.absoluteMomentum();
                if (p2 < p1/2) {
                    scatteredParameters[Acts::eBoundQOverP] =
                            beforeParticle.charge()/(0.999*beforeParticle.absoluteMomentum());
                } else {
                    scatteredParameters[Acts::eBoundQOverP] =
                            tempParticle.charge()/tempParticle.absoluteMomentum();
                }
            }

            Acts::SquareMatrix2 cov = Acts::SquareMatrix2::Identity();
            Acts::BoundSquareMatrix resetCov = Acts::BoundSquareMatrix::Identity();
            scatteredParameters[Acts::eBoundPhi] = tempParticle.phi();
            scatteredParameters[Acts::eBoundTheta] = tempParticle.theta();

            stepper.resetState(state.stepping, scatteredParameters, resetCov, surface);
            SimpleSourceLink sl(loc, cov, geoId, sourceId);
            LUXEDataContainer::SimMeasurement SM{sl, scatteredParameters, trueVertex, static_cast<int>(sourceId)};
            result.push_back(SM);

        }
    };
/// Propagate the track create smeared measurements from local coordinates.

    class MeasurementsCreator : public IAlgorithm {
    public:

        /// @brief The nested configuration struct
        struct Config {
            /// The output collection
            std::string outputCollection = "Measurements";
            // detector
            std::shared_ptr<const Acts::Experimental::Detector> detector;
            // variable bins magnetic field
            std::shared_ptr<Acts::InterpolatedBFieldMap<LUXEMagneticField::vGrid>> BFieldPtr;
        };

        /// @brief Constructor
        MeasurementsCreator(Config config, Acts::Logging::Level level)
                : IAlgorithm("MeasurementsCreator", level),
                  m_cfg(std::move(config)) {
            m_outputMeasurements.initialize(m_cfg.outputCollection);
        }
        ~MeasurementsCreator() = default;

        template <typename propagator_t, typename track_parameters_t>
        LUXEDataContainer::SimMeasurements createMeasurements(const propagator_t& propagator,
                                       const AlgorithmContext& ctx,
                                       const track_parameters_t& trackParameters) const {
            using namespace Acts::UnitLiterals;
            using Actions = Acts::ActionList<LUXENavigator::MeasurementsCreatorAction>;
            using Aborters = Acts::AbortList<Acts::EndOfWorldReached>;
            // Set options for propagator
            Acts::PropagatorOptions<Actions, Aborters> options(ctx.geoContext, ctx.magFieldContext);

            auto& creator = options.actionList.get<LUXENavigator::MeasurementsCreatorAction>();
            creator.sourceId = ctx.eventNumber;

            // Launch and collect the measurements
            auto result = propagator.propagate(trackParameters, options).value();
            return std::move(result.template get<LUXEDataContainer::SimMeasurements>());
        }

        /// @brief The execute method
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            using namespace Acts::UnitLiterals;

            // Create IP covariance matrix from
            // reasonable standard deviations
            Acts::BoundVector ipStdDev;
            ipStdDev[Acts::eBoundLoc0] = 100_um;
            ipStdDev[Acts::eBoundLoc1] = 100_um;
            ipStdDev[Acts::eBoundTime] = 25_ns;
            ipStdDev[Acts::eBoundPhi] = 2_degree;
            ipStdDev[Acts::eBoundTheta] = 2_degree;
            ipStdDev[Acts::eBoundQOverP] = 1 / 100_GeV;
            Acts::BoundSquareMatrix ipCov =
                    ipStdDev.cwiseProduct(ipStdDev).asDiagonal();

            // Create the measurements
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<> pDisP(0.002,0.0018);
            std::normal_distribution<> pDisM(-0.002,0.0018);
            std::gamma_distribution<double> pzDis(3, 1.2);
            std::uniform_real_distribution<> uni(2.2,2.3);

            Acts::ActsScalar m_e = 0.000511;
            LUXEDataContainer::SimMeasurements results;
            Acts::ActsScalar px = (pDisP(gen)+pDisM(gen))/2;
            Acts::ActsScalar pz = (pDisP(gen)+pDisM(gen))/2;
            Acts::ActsScalar py = pzDis(gen)+1;
//              Acts::ActsScalar py = uni(gen);
            Acts::ActsScalar p = std::sqrt(std::pow(px,2)+std::pow(py,2)+std::pow(pz,2));
            Acts::ActsScalar E = std::hypot(p,m_e);
            Acts::ActsScalar theta = std::acos(pz / p);
            Acts::ActsScalar phi = std::atan2(py, px);
            Acts::Vector4 mPos4{0,0,0,0};

            Acts::Experimental::DetectorNavigator::Config cfg{&(*m_cfg.detector)};
            cfg.resolvePassive = false;
            cfg.resolveMaterial = true;
            cfg.resolveSensitive = true;
            Acts::Experimental::DetectorNavigator navigator(
                    cfg, Acts::getDefaultLogger("Detector Navigation", Acts::Logging::VERBOSE));
            Acts::EigenStepper<> stepper(std::move(m_cfg.BFieldPtr));

            auto propagator = Acts::Propagator<decltype(stepper), Acts::Experimental::DetectorNavigator>(
                    std::move(stepper), std::move(navigator));
            results = createMeasurements(
                    propagator, ctx,
                    Acts::CurvilinearTrackParameters(mPos4, phi*180/M_PI*1_degree, theta*180/M_PI*1_degree,
                                                     1_e / (p*1_GeV), ipCov, Acts::ParticleHypothesis::electron()));

            m_outputMeasurements(ctx, std::move(results));

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;

        WriteDataHandle<LUXEDataContainer::SimMeasurements> m_outputMeasurements
                {this, "OutputMeasurements"};
    };
} // LUXENavigator