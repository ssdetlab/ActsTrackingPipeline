#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/EigenStepperError.hpp"
#include "Acts/Propagator/PropagatorTraits.hpp"
#include "Acts/Propagator/StepperExtensionList.hpp"
#include "Acts/Utilities/QuickMath.hpp"
#include "Acts/Utilities/Result.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

#include "TrackingPipeline/Propagation/FieldGradientExtension.hpp"

/// @brief Runge-Kutta-Nystroem stepper with the magnetic field gradient
/// in the covariance transport
///
/// Acts::EigenStepper with the FieldGradientExtension. The step is
/// Acts::EigenStepper::step, except for the accumulation of the
/// transport jacobian, which uses the full matrix product.
class FieldGradientEigenStepper
    : public Acts::EigenStepper<
          Acts::StepperExtensionList<FieldGradientExtension>> {
 public:
  /// Base stepper shorthand
  using Base =
      Acts::EigenStepper<Acts::StepperExtensionList<FieldGradientExtension>>;

  /// Constructors of the base stepper
  using Base::Base;

  /// @brief Get the field gradient for the stepping
  ///
  /// The derivative is left unchanged by the field providers
  /// that do not evaluate it.
  ///
  /// @param state stepping state
  /// @param pos global position of the field inquiry
  /// @param derivative field jacobian in global coordinates,
  /// (i, j) = dB_i / dx_j
  ///
  /// @return field value in the result wrapper
  Acts::Result<Acts::Vector3> getFieldGradient(
      State& state, const Acts::Vector3& pos,
      Acts::ActsMatrix<3, 3>& derivative) const {
    return m_bField->getFieldGradient(pos, derivative, state.fieldCache);
  }

  /// @brief Perform a Runge-Kutta track parameter propagation step
  ///
  /// Hides Acts::EigenStepper::step.
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param state propagation state
  /// @param navigator navigator of the propagation
  ///
  /// @return step size in the result wrapper
  template <typename propagator_state_t, typename navigator_t>
  Acts::Result<double> step(propagator_state_t& state,
                            const navigator_t& navigator) const {
    // Adapted from Acts::EigenStepper::step (Acts v36,
    // Core/include/Acts/Propagator/EigenStepper.ipp),
    // Copyright (C) 2019-2022 CERN for the benefit of the Acts project,
    // Mozilla Public License, v. 2.0, http://mozilla.org/MPL/2.0/.

    // Runge-Kutta integrator state
    auto& sd = state.stepping.stepData;

    double errorEstimate = 0;
    double h2 = 0;
    double half_h = 0;

    auto pos = position(state.stepping);
    auto dir = direction(state.stepping);

    // First Runge-Kutta point (at current position)
    auto fieldRes = getField(state.stepping, pos);
    if (!fieldRes.ok()) {
      return fieldRes.error();
    }
    sd.B_first = *fieldRes;
    if (!state.stepping.extension.validExtensionForStep(state, *this,
                                                        navigator) ||
        !state.stepping.extension.k1(state, *this, navigator, sd.k1, sd.B_first,
                                     sd.kQoP)) {
      return 0.;
    }

    const auto calcStepSizeScaling =
        [&](const double errorEstimate_) -> double {
      // For details about these values see ATL-SOFT-PUB-2009-001
      constexpr double lower = 0.25;
      constexpr double upper = 4.0;
      // This is given by the order of the Runge-Kutta method
      constexpr double exponent = 0.25;

      // Whether to use fast power function if available
      constexpr bool tryUseFastPow{false};

      double x = state.options.stepping.stepTolerance / errorEstimate_;

      if constexpr (exponent == 0.25 && !tryUseFastPow) {
        // This is 3x faster than std::pow
        x = std::sqrt(std::sqrt(x));
      } else if constexpr (std::numeric_limits<double>::is_iec559 &&
                           tryUseFastPow) {
        x = Acts::fastPow(x, exponent);
      } else {
        x = std::pow(x, exponent);
      }

      return std::clamp(x, lower, upper);
    };

    const auto isErrorTolerable = [&](const double errorEstimate_) {
      // For details about these values see ATL-SOFT-PUB-2009-001
      constexpr double marginFactor = 4.0;

      return errorEstimate_ <=
             marginFactor * state.options.stepping.stepTolerance;
    };

    // The following functor starts to perform a Runge-Kutta step of a certain
    // size, going up to the point where it can return an estimate of the local
    // integration error. The results are stated in the local variables above,
    // allowing integration to continue once the error is deemed satisfactory
    const auto tryRungeKuttaStep = [&](const double h) -> Acts::Result<bool> {
      // helpers because bool and std::error_code are ambiguous
      constexpr auto success = &Acts::Result<bool>::success;
      constexpr auto failure = &Acts::Result<bool>::failure;

      // State the square and half of the step size
      h2 = h * h;
      half_h = h * 0.5;

      // Second Runge-Kutta point
      const Acts::Vector3 pos1 = pos + half_h * dir + h2 * 0.125 * sd.k1;
      auto field = getField(state.stepping, pos1);
      if (!field.ok()) {
        return failure(field.error());
      }
      sd.B_middle = *field;

      if (!state.stepping.extension.k2(state, *this, navigator, sd.k2,
                                       sd.B_middle, sd.kQoP, half_h, sd.k1)) {
        return success(false);
      }

      // Third Runge-Kutta point
      if (!state.stepping.extension.k3(state, *this, navigator, sd.k3,
                                       sd.B_middle, sd.kQoP, half_h, sd.k2)) {
        return success(false);
      }

      // Last Runge-Kutta point
      const Acts::Vector3 pos2 = pos + h * dir + h2 * 0.5 * sd.k3;
      field = getField(state.stepping, pos2);
      if (!field.ok()) {
        return failure(field.error());
      }
      sd.B_last = *field;
      if (!state.stepping.extension.k4(state, *this, navigator, sd.k4,
                                       sd.B_last, sd.kQoP, h, sd.k3)) {
        return success(false);
      }

      // Compute and check the local integration error estimate
      errorEstimate =
          h2 * ((sd.k1 - sd.k2 - sd.k3 + sd.k4).template lpNorm<1>() +
                std::abs(sd.kQoP[0] - sd.kQoP[1] - sd.kQoP[2] + sd.kQoP[3]));
      // Protect against division by zero
      errorEstimate = std::max(1e-20, errorEstimate);

      return success(isErrorTolerable(errorEstimate));
    };

    const double initialH =
        state.stepping.stepSize.value() * state.options.direction;
    double h = initialH;
    std::size_t nStepTrials = 0;
    // Select and adjust the appropriate Runge-Kutta step size as given
    // ATL-SOFT-PUB-2009-001
    while (true) {
      nStepTrials++;
      auto res = tryRungeKuttaStep(h);
      if (!res.ok()) {
        return res.error();
      }
      if (!!res.value()) {
        break;
      }

      const double stepSizeScaling = calcStepSizeScaling(errorEstimate);
      h *= stepSizeScaling;

      // If step size becomes too small the particle remains at the initial
      // place
      if (std::abs(h) < std::abs(state.options.stepping.stepSizeCutOff)) {
        // Not moving due to too low momentum needs an aborter
        return Acts::EigenStepperError::StepSizeStalled;
      }

      // If the parameter is off track too much or given stepSize is not
      // appropriate
      if (nStepTrials > state.options.stepping.maxRungeKuttaStepTrials) {
        // Too many trials, have to abort
        return Acts::EigenStepperError::StepSizeAdjustmentFailed;
      }
    }

    // When doing error propagation, update the associated Jacobian matrix
    if (state.stepping.covTransport) {
      // using the direction before updated below

      // The step transport matrix in global coordinates
      Acts::FreeMatrix D;
      if (!state.stepping.extension.finalize(state, *this, navigator, h, D)) {
        return Acts::EigenStepperError::StepInvalid;
      }

      // Acts::EigenStepper accumulates D with a blocked product that
      // assumes the top left 4x4 block of D to be identity and the bottom
      // left one to be zero (dF/dx = I, dG/dx = 0). With the field gradient
      // both blocks are non-trivial, so the full product is used.
      state.stepping.jacTransport = D * state.stepping.jacTransport;
    } else {
      if (!state.stepping.extension.finalize(state, *this, navigator, h)) {
        return Acts::EigenStepperError::StepInvalid;
      }
    }

    // Update the track parameters according to the equations of motion
    state.stepping.pars.template segment<3>(Acts::eFreePos0) +=
        h * dir + h2 / 6. * (sd.k1 + sd.k2 + sd.k3);
    state.stepping.pars.template segment<3>(Acts::eFreeDir0) +=
        h / 6. * (sd.k1 + 2. * (sd.k2 + sd.k3) + sd.k4);
    (state.stepping.pars.template segment<3>(Acts::eFreeDir0)).normalize();

    if (state.stepping.covTransport) {
      // using the updated direction
      state.stepping.derivative.template head<3>() =
          state.stepping.pars.template segment<3>(Acts::eFreeDir0);
      state.stepping.derivative.template segment<3>(4) = sd.k4;
    }

    state.stepping.pathAccumulated += h;
    ++state.stepping.nSteps;
    state.stepping.nStepTrials += nStepTrials;

    const double stepSizeScaling = calcStepSizeScaling(errorEstimate);
    const double nextAccuracy = std::abs(h * stepSizeScaling);
    const double previousAccuracy =
        std::abs(state.stepping.stepSize.accuracy());
    const double initialStepLength = std::abs(initialH);
    if (nextAccuracy < initialStepLength || nextAccuracy > previousAccuracy) {
      state.stepping.stepSize.setAccuracy(nextAccuracy);
    }

    return h;
  }
};

namespace Acts {

/// The stepper supports the bound track parameters
template <typename navigator_t>
struct SupportsBoundParameters<FieldGradientEigenStepper, navigator_t>
    : public std::true_type {};

}  // namespace Acts
