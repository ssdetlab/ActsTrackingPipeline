#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <random>
#include <vector>

/// Propagator action to create measurements.
namespace LUXENavigator {
enum class MeasurementType {
    eLoc0,
    eLoc1,
    eLoc01,
};
struct MeasurementResolution {
    MeasurementType type = MeasurementType::eLoc0;
    // Depending on the type, only the first value is used.
    std::array<double, 2> stddev = {50 * Acts::UnitConstants::um, 0};
};

/// Measurement resolution configuration for a full detector geometry.
using MeasurementResolutionMap =
        Acts::GeometryHierarchyMap<MeasurementResolution>;

/// Result struct for generated measurements and outliers.
struct Measurements {
    std::vector<Acts::detail::Test::TestSourceLink> sourceLinks;
    std::vector<Acts::BoundVector> truthParameters;
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
        std::cout<<"check4"<<std::endl;
        // Apply global to local
        Acts::Vector2 loc =
                surface.globalToLocal(state.geoContext, stepper.position(state.stepping),
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
        result.truthParameters.push_back(std::move(parameters));

        Acts::SquareMatrix2 cov = Acts::SquareMatrix2::Identity();

        Acts::Vector2 val = loc;

        result.sourceLinks.emplace_back(Acts::eBoundLoc0, Acts::eBoundLoc1, val, cov, geoId,
                                        sourceId);
        result.truthParameters.push_back(parameters);
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
    options.pathLimit = 6_m;

    // Launch and collect the measurements
    auto result = propagator.propagate(trackParameters, options).value();
    return std::move(result.template get<Measurements>());
}
} // LUXENavigator