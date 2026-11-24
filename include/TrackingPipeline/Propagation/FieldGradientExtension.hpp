#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Propagator/DefaultExtension.hpp"

#include <cmath>
#include <cstddef>

/// @brief RKN4 stepper extension with the magnetic field gradient
/// in the transport matrix
///
/// The k_i's are evaluated by Acts::DefaultExtension. The transport
/// matrix D of the step additionally contains the dependence of the
/// k_i's on the position through the field gradient dB/dx, evaluated
/// at the three field evaluation points of the step. With it, dF/dx
/// and dG/dx are not identity and zero, as Acts::DefaultExtension
/// assumes, so D has to be accumulated with the full matrix product
/// (see FieldGradientEigenStepper).
///
/// The stepper has to provide
/// getFieldGradient(state, position, derivative). Field providers
/// that do not fill the derivative contribute a zero gradient.
struct FieldGradientExtension : public Acts::DefaultExtension {
  /// @brief Veto function after a RKN4 step was accepted by judging on the
  /// error of the step. Propagates the time.
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param state State of the propagator
  /// @param stepper Stepper of the propagation
  /// @param navigator Navigator of the propagation
  /// @param h Step size
  ///
  /// @return Boolean flag if the calculation is valid
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool finalize(propagator_state_t& state, const stepper_t& stepper,
                const navigator_t& /*navigator*/, const double h) const {
    propagateTime(state, stepper, h);
    return true;
  }

  /// @brief Veto function after a RKN4 step was accepted by judging on the
  /// error of the step. Propagates the time and evaluates the transport
  /// matrix.
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param state State of the propagator
  /// @param stepper Stepper of the propagation
  /// @param navigator Navigator of the propagation
  /// @param h Step size
  /// @param D Transport matrix
  ///
  /// @return Boolean flag if the calculation is valid
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool finalize(propagator_state_t& state, const stepper_t& stepper,
                const navigator_t& /*navigator*/, const double h,
                Acts::FreeMatrix& D) const {
    propagateTime(state, stepper, h);
    return transportMatrix(state, stepper, h, D);
  }

 private:
  /// Free parameters jacobian of a 3-vector
  using FreeJacobian3 = Acts::ActsMatrix<3, Acts::eFreeSize>;

  /// @brief Propagation function for the time coordinate,
  /// dt/ds = sqrt(1 + m^2 / p^2)
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param state State of the propagator
  /// @param stepper Stepper of the propagation
  /// @param h Step size
  template <typename propagator_state_t, typename stepper_t>
  void propagateTime(propagator_state_t& state, const stepper_t& stepper,
                     const double h) const {
    auto m = stepper.particleHypothesis(state.stepping).mass();
    auto p = stepper.absoluteMomentum(state.stepping);
    auto dtds = std::sqrt(1 + m * m / (p * p));
    state.stepping.pars[Acts::eFreeTime] += h * dtds;
    if (state.stepping.covTransport) {
      state.stepping.derivative(3) = dtds;
    }
  }

  /// @brief Derivative of qop * (a x B(x)) with respect to the free
  /// parameters at the start of the step,
  /// qop * (da x B + a x (dB/dx * dx)) + dqop * (a x B)
  ///
  /// @param qop charge over momentum
  /// @param a direction argument
  /// @param da derivative of the direction argument
  /// @param bField field value at the evaluation point
  /// @param dBdx field gradient at the evaluation point
  /// @param dx derivative of the evaluation point
  ///
  /// @return derivative of the k_i
  static FreeJacobian3 kDerivative(double qop, const Acts::Vector3& a,
                                   const FreeJacobian3& da,
                                   const Acts::Vector3& bField,
                                   const Acts::ActsMatrix<3, 3>& dBdx,
                                   const FreeJacobian3& dx) {
    FreeJacobian3 dk;
    for (std::size_t i = 0; i < Acts::eFreeSize; i++) {
      Acts::Vector3 dai = da.col(i);
      Acts::Vector3 dxi = dx.col(i);
      dk.col(i) = qop * (dai.cross(bField) + a.cross(dBdx * dxi));
    }
    dk.col(Acts::eFreeQOverP) += a.cross(bField);
    return dk;
  }

  /// @brief Calculates the transport matrix D of the step,
  /// differentiating the RKN4 step of Acts::EigenStepper
  /// with the field gradient
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param state State of the propagator
  /// @param stepper Stepper of the propagation
  /// @param h Step size
  /// @param D Transport matrix
  ///
  /// @return Boolean flag if evaluation is valid
  template <typename propagator_state_t, typename stepper_t>
  bool transportMatrix(propagator_state_t& state, const stepper_t& stepper,
                       const double h, Acts::FreeMatrix& D) const {
    auto m = state.stepping.particleHypothesis.mass();
    auto& sd = state.stepping.stepData;
    Acts::Vector3 pos = stepper.position(state.stepping);
    Acts::Vector3 dir = stepper.direction(state.stepping);
    double qop = stepper.qOverP(state.stepping);
    auto p = stepper.absoluteMomentum(state.stepping);
    auto dtds = std::sqrt(1 + m * m / (p * p));

    double half_h = h * 0.5;
    double h2 = h * h;

    // Field evaluation points of the step, as in Acts::EigenStepper::step
    Acts::Vector3 posMiddle = pos + half_h * dir + h2 * 0.125 * sd.k1;
    Acts::Vector3 posLast = pos + h * dir + h2 * 0.5 * sd.k3;

    // Field gradients at the evaluation points
    Acts::ActsMatrix<3, 3> dBdxFirst = Acts::ActsMatrix<3, 3>::Zero();
    Acts::ActsMatrix<3, 3> dBdxMiddle = Acts::ActsMatrix<3, 3>::Zero();
    Acts::ActsMatrix<3, 3> dBdxLast = Acts::ActsMatrix<3, 3>::Zero();
    if (!stepper.getFieldGradient(state.stepping, pos, dBdxFirst).ok() ||
        !stepper.getFieldGradient(state.stepping, posMiddle, dBdxMiddle).ok() ||
        !stepper.getFieldGradient(state.stepping, posLast, dBdxLast).ok()) {
      return false;
    }

    // Derivatives of the start position and direction
    FreeJacobian3 dPos = FreeJacobian3::Zero();
    dPos.block<3, 3>(0, Acts::eFreePos0) = Acts::ActsMatrix<3, 3>::Identity();
    FreeJacobian3 dDir = FreeJacobian3::Zero();
    dDir.block<3, 3>(0, Acts::eFreeDir0) = Acts::ActsMatrix<3, 3>::Identity();

    // Derivatives of the k_i's and of the field evaluation points
    FreeJacobian3 dk1 =
        kDerivative(qop, dir, dDir, sd.B_first, dBdxFirst, dPos);

    FreeJacobian3 dPosMiddle = dPos + half_h * dDir + h2 * 0.125 * dk1;
    FreeJacobian3 dk2 =
        kDerivative(qop, dir + half_h * sd.k1, dDir + half_h * dk1, sd.B_middle,
                    dBdxMiddle, dPosMiddle);
    FreeJacobian3 dk3 =
        kDerivative(qop, dir + half_h * sd.k2, dDir + half_h * dk2, sd.B_middle,
                    dBdxMiddle, dPosMiddle);

    FreeJacobian3 dPosLast = dPos + h * dDir + h2 * 0.5 * dk3;
    FreeJacobian3 dk4 = kDerivative(qop, dir + h * sd.k3, dDir + h * dk3,
                                    sd.B_last, dBdxLast, dPosLast);

    D = Acts::FreeMatrix::Identity();

    // dF/du, with F the position update
    D.block<3, Acts::eFreeSize>(Acts::eFreePos0, 0) =
        dPos + h * dDir + h2 / 6. * (dk1 + dk2 + dk3);

    // dG/du, with G the direction update
    D.block<3, Acts::eFreeSize>(Acts::eFreeDir0, 0) =
        dDir + h / 6. * (dk1 + 2. * (dk2 + dk3) + dk4);

    D(Acts::eFreeTime, Acts::eFreeQOverP) = h * m * m * qop / dtds;
    return true;
  }
};
