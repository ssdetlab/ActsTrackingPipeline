#pragma once

#include "Acts/Definitions/Units.hpp"

#include <cmath>

using namespace Acts::UnitLiterals;

namespace detail {

/// @brief calculate Theta RMS of multiple scattering
///
/// @param X0 radiation length of the material [g / cm2]
/// @param rho densoty of the material [g / cm3]
/// @param x thickness of the material [cm]
/// @param P particle momentum [MeV]
/// @param z particle charge [e]
///
/// @return RMS of the multiple scattering angle
double getMcpThetaRms(double X0, double rho, double x, double P, double z) {
  double t = rho * x;
  double thetaRms = 13.6_MeV / P * z * std::sqrt(t / X0) *
                    (1 + 0.038 * std::log(t * z * z / (X0)));
  return thetaRms;
}

/// @brief specialization for Silicon
///
/// @param x thickness of the material [cm]
/// @param P particle momentum [MeV]
/// @param z particle charge [e]
///
/// @return RMS of the multiple scattering angle in Si
double getMcpThetaRmsSi(double x, double P, double z) {
  double X0 = 21.82 * Acts::UnitConstants::g / Acts::UnitConstants::cm2;
  double rho = 2.329 * Acts::UnitConstants::g / Acts::UnitConstants::cm3;
  return getMcpThetaRms(X0, rho, x, P, z);
}

}  // namespace detail
