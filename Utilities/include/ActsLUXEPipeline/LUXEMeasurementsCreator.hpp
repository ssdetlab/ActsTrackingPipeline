#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/Scattering.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/BetheBloch.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/BetheHeitler.hpp"

#include "Acts/Material/Material.hpp"


#include <memory>
#include <random>
#include <vector>
#include <iostream>

/// Propagator action to create measurements.
namespace LUXENavigator {

// utility function to build a particle from the dataset parameters
ActsFatras::Particle makeParticle(Acts::PdgParticle pdg, Acts::Vector3 dir, double p, Acts::Vector4 pos4) {
    const auto id = ActsFatras::Barcode().setVertexPrimary(1).setParticle(1);
    return ActsFatras::Particle(id, pdg)
            .setPosition4(pos4)
            .setDirection(dir)
            .setAbsoluteMomentum(p);
}
//    .setDirection(std::sin(theta) * std::cos(phi),
//            std::sin(theta) * std::sin(phi), std::cos(theta))
Acts::MaterialSlab makeSiliconSlab() {
    using namespace Acts::UnitLiterals;
    return {Acts::Material::fromMolarDensity(9.370_cm, 46.52_cm, 28.0855, 14,
                                      (2.329 / 28.0855) * 1_mol / 1_cm3), 100_um};
}

enum class MeasurementType {
    eLoc0,
    eLoc1,
    eLoc01,
};
struct MeasurementResolution {
    MeasurementType type = MeasurementType::eLoc0;
    // Depending on the type, only the first value is used.
    std::array<double, 2> stddev = {15 * Acts::UnitConstants::um, 15 * Acts::UnitConstants::um};
};

/// Measurement resolution configuration for a full detector geometry.
using MeasurementResolutionMap =
        Acts::GeometryHierarchyMap<MeasurementResolution>;

/// Result struct for generated measurements and outliers.
struct Measurements {
    unsigned int eventId;
    std::vector<SimpleSourceLink> sourceLinks;
    std::vector<Acts::Vector3> fullTrack;
    std::vector<Acts::BoundVector> truthParameters;
    std::vector<Acts::Vector3> globalPosition;
};

struct MeasurementsCreator {
    using result_type = Measurements;

    MeasurementResolutionMap resolutions;
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

        // only generate measurements on surfaces
        if (!navigator.currentSurface(state.navigation)) {
            result.fullTrack.push_back(stepper.position(state.stepping));
            return;
        }
        const Acts::Surface &surface = *navigator.currentSurface(state.navigation);
        const Acts::GeometryIdentifier geoId = surface.geometryId();
        // only generate measurements on sensitive surface
        if (!geoId.sensitive()) {
            ACTS_VERBOSE("Create no measurements on non-sensitive surface " << geoId);
            return;
        }
        // only generate measurements if a resolution is configured
        auto found = resolutions.find(geoId);
        if (found == resolutions.end()) {
            ACTS_VERBOSE("No resolution configured for sensitive surface " << geoId);
            return;
        }
        const MeasurementResolution &resolution = *found;
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

        result.eventId = sourceId;
        result.sourceLinks.emplace_back(loc, cov, geoId, sourceId);
        result.truthParameters.push_back(std::move(scatteredParameters));
        result.globalPosition.push_back(std::move(pos3));
    }
};
/// Propagate the track create smeared measurements from local coordinates.
template <typename propagator_t, typename track_parameters_t>
Measurements createMeasurements(const propagator_t& propagator,
                                const Acts::GeometryContext& geoCtx,
                                const Acts::MagneticFieldContext& magCtx,
                                const track_parameters_t& trackParameters,
                                const MeasurementResolutionMap& resolutions,
                                std::size_t sourceId = 0u) {
    using namespace Acts::UnitLiterals;
    using Actions = Acts::ActionList<LUXENavigator::MeasurementsCreator>;
    using Aborters = Acts::AbortList<Acts::EndOfWorldReached>;

    // Set options for propagator
    Acts::PropagatorOptions<Actions, Aborters> options(geoCtx, magCtx);

    auto& creator = options.actionList.get<LUXENavigator::MeasurementsCreator>();
    creator.resolutions = resolutions;
    creator.sourceId = sourceId;

    // Launch and collect the measurements
    auto result = propagator.propagate(trackParameters, options).value();
    return std::move(result.template get<Measurements>());
}

} // LUXENavigator

// COVARIANCE AVE

//stepper.transportCovarianceToBound(state.stepping, surface,
//Acts::FreeToBoundCorrection(true));
//std::random_device rd;
//std::mt19937 rng(rd());
//std::normal_distribution<double> normalDist(0., 1.);
//
//Acts::Vector2 stddev(resolution.stddev[0], resolution.stddev[1]);
//Acts::SquareMatrix2 cov = stddev.cwiseProduct(stddev).asDiagonal();
//std::cout<<state.stepping.cov<<std::endl;
//Acts::Vector2 val = loc + stddev.cwiseProduct(
//        Acts::Vector2(normalDist(rng), normalDist(rng)));