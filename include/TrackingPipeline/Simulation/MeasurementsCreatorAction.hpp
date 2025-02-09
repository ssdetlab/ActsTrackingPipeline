#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/BetheBloch.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/BetheHeitler.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/Scattering.hpp"
#include <Acts/Definitions/Common.hpp>

#include <random>
#include <stdexcept>
#include <vector>

using namespace Acts::UnitLiterals;

/// @brief Action type that creates measurements on sensitive surfaces
///
/// Action type that is called by the navigator to create measurements
/// on the sensitive surfaces, when encountered. Takes into account
/// multiple scattering and energy loss processes, if the surfaces
/// has material assigned
struct MeasurementsCreatorAction {
  using result_type = std::vector<Acts::BoundTrackParameters>;

  std::int32_t sourceId = 0;

  /// @brief Operator that is callable by an ActionList. The function
  /// collects the surfaces
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [out] result Vector of matching surfaces
  /// @param [in] state State of the propagator
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  void operator()(propagator_state_t &state, const stepper_t &stepper,
                  const navigator_t &navigator, result_type &result,
                  const Acts::Logger &logger) const {
    // Only generate measurements on surfaces
    if (!navigator.currentSurface(state.navigation)) {
      return;
    }
    const Acts::Surface *surface = navigator.currentSurface(state.navigation);
    const Acts::GeometryIdentifier &geoId = surface->geometryId();

    // only generate measurements on sensitive surface
    if (!geoId.sensitive()) {
      ACTS_VERBOSE("Create no measurements on non-sensitive surface " << geoId);
      return;
    }

    // Apply global to local
    Acts::Vector3 globalPos = stepper.position(state.stepping);
    Acts::Vector2 localPos =
        surface
            ->globalToLocal(state.geoContext, globalPos,
                            stepper.direction(state.stepping))
            .value();

    // The truth info of local hit
    Acts::BoundVector parameters = Acts::BoundVector::Zero();
    parameters[Acts::eBoundLoc0] = localPos[Acts::eBoundLoc0];
    parameters[Acts::eBoundLoc1] = localPos[Acts::eBoundLoc1];

    const Acts::Vector3 &direction = stepper.direction(state.stepping);
    parameters[Acts::eBoundPhi] = Acts::VectorHelpers::phi(direction);
    parameters[Acts::eBoundTheta] = Acts::VectorHelpers::theta(direction);
    parameters[Acts::eBoundQOverP] = state.stepping.pars[Acts::eFreeQOverP];
    parameters[Acts::eBoundTime] = state.stepping.pars[Acts::eFreeTime];

    // Construct a particle object
    Acts::Vector4 globalFourPos = {globalPos.x(), globalPos.y(), globalPos.z(),
                                   parameters[Acts::eBoundTime]};

    const auto id = ActsFatras::Barcode().setVertexPrimary(1).setParticle(1);
    auto simParticle = ActsFatras::Particle(
        id, stepper.particleHypothesis(state.stepping).absolutePdg());
    simParticle.setPosition4(globalFourPos);
    simParticle.setDirection(direction);
    simParticle.setAbsoluteMomentum(stepper.charge(state.stepping) /
                                    parameters[Acts::eBoundQOverP]);
    simParticle.setReferenceSurface(surface);

    // Create the scattering and energy loss processes
    std::random_device rd;
    std::mt19937 gen(rd());

    // Scatter particle parameters
    Acts::BoundVector scatteredParameters = parameters;

    if (surface->surfaceMaterial() != nullptr) {
      // Retrieve the material
      Acts::MaterialSlab material;
      auto bsm = dynamic_cast<const Acts::BinnedSurfaceMaterial *>(
          surface->surfaceMaterial());
      auto hsm = dynamic_cast<const Acts::HomogeneousSurfaceMaterial *>(
          surface->surfaceMaterial());
      auto psm = dynamic_cast<const Acts::ProtoSurfaceMaterial *>(
          surface->surfaceMaterial());

      if (bsm != nullptr) {
        material = bsm->materialSlab(localPos);
      } else if (hsm != nullptr) {
        material = hsm->materialSlab(localPos);
      } else if (psm != nullptr) {
        material = psm->materialSlab(localPos);
      } else {
        throw std::runtime_error("Unsupported surface material encountered");
      }

      auto scattering = ActsFatras::GaussianMixtureScattering();
      ActsFatras::BetheBloch BBProcess;
      ActsFatras::BetheHeitler BHProcess;

      // Apply the scattering and energy loss processes
      auto scatter = scattering(gen, material, simParticle);
      auto BBLoss = BBProcess(gen, material, simParticle);
      auto BHLoss = BHProcess(gen, material, simParticle);
      scatteredParameters[Acts::eBoundQOverP] =
          simParticle.charge() / simParticle.absoluteMomentum();

      scatteredParameters[Acts::eBoundPhi] = simParticle.phi();
      scatteredParameters[Acts::eBoundTheta] =
          (std::abs(std::sin(simParticle.theta())) > 1e-6) ? simParticle.theta()
                                                           : 1e-6;
    }

    // Reset the state
    stepper.resetState(state.stepping, scatteredParameters, state.stepping.cov,
                       *surface);

    Acts::BoundTrackParameters boundsPars(
        surface->getSharedPtr(), scatteredParameters, state.stepping.cov,
        stepper.particleHypothesis(state.stepping));
    result.push_back(boundsPars);
  }
};
