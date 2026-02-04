#include "TrackingPipeline/Simulation/SphericalMomentumGenerator.hpp"

#include "Acts/Definitions/Algebra.hpp"

SphericalMomentumGenerator::SphericalMomentumGenerator(const Config& config)
    : m_cfg(config) {
  double pMin = m_cfg.pRange.first;
  double pMax = m_cfg.pRange.second;
  double meanP = (pMax + pMin) / 2.0;
  double varP = (pMax - pMin) * (pMax - pMin) / 12.0;

  bool nonBijective = false;

  double sinThetaMin = std::sin(m_cfg.thetaRange.first);
  double cosThetaMin = std::cos(m_cfg.thetaRange.first);

  double sinThetaMax = std::sin(m_cfg.thetaRange.second);
  double cosThetaMax = std::cos(m_cfg.thetaRange.second);

  double rhoMin = std::min(sinThetaMin, sinThetaMax);
  double rhoMax = std::max(sinThetaMin, sinThetaMax);

  double rhoMin2 = rhoMin * rhoMin;
  double rhoMax2 = rhoMax * rhoMax;
  if (m_cfg.thetaRange.first < M_PI_2 && m_cfg.thetaRange.second > M_PI_2) {
    nonBijective = true;
  }

  double zMin = std::min(cosThetaMin, cosThetaMax);
  double zMax = std::max(cosThetaMin, cosThetaMax);
  double meanZ = (zMax + zMin) / 2.0;
  double varZ = (zMax - zMin) * (zMax - zMin) / 12.0;

  double phiMin = m_cfg.phiRange.first;
  double phiMax = m_cfg.phiRange.second;
  double meanPhi = (phiMax + phiMin) / 2.0;
  double varPhi = (phiMax - phiMin) * (phiMax - phiMin) / 12.0;

  double meanRho;
  double mean2Rho;
  if (!nonBijective) {
    meanRho =
        1.0 / 2.0 *
        (rhoMin * std::sqrt(1 - rhoMin2) - rhoMax * std::sqrt(1 - rhoMax2) -
         std::asin(rhoMin) + std::asin(rhoMax)) /
        (std::sqrt(1 - rhoMin2) - std::sqrt(1 - rhoMax2));

    mean2Rho =
        1.0 / 3.0 *
        (1 + rhoMin2 + rhoMax2 - std::sqrt((-1 + rhoMin2) * (-1 + rhoMax2)));
  } else {
    double weight1 = std::max(zMax, std::abs(zMin)) / (zMax - zMin);
    double weight2 = std::min(zMax, std::abs(zMin)) / (zMax - zMin);

    double meanRho1 =
        1.0 / 2.0 *
        (rhoMin * std::sqrt(1 - rhoMin2) - std::asin(rhoMin) + M_PI_2) /
        std::sqrt(1 - rhoMin2);

    double meanRho2 =
        1.0 / 2.0 *
        (rhoMax * std::sqrt(1 - rhoMax2) - std::asin(rhoMax) + M_PI_2) /
        std::sqrt(1 - rhoMax2);

    double mean2Rho1 = 1.0 / 3.0 * (2 + rhoMin2);

    double mean2Rho2 = 1.0 / 3.0 * (2 + rhoMax2);

    meanRho = meanRho1 * weight1 + meanRho2 * weight2;
    mean2Rho = mean2Rho1 * weight1 + mean2Rho2 * weight2;
  }
  double varR = mean2Rho - meanRho * meanRho;

  double meanZR = (-std::pow(1 - zMax * zMax, 3.0 / 2.0) +
                   std::pow(1 - zMin * zMin, 3.0 / 2.0)) /
                  (3.0 * (zMax - zMin));
  double covZR = meanZR - meanZ * meanRho;

  Acts::SquareMatrix3 covCyl = Acts::SquareMatrix3::Zero();
  covCyl(0, 0) = varR;
  covCyl(0, 2) = covZR;

  covCyl(1, 1) = varPhi;

  covCyl(2, 0) = covZR;
  covCyl(2, 2) = varZ;

  Acts::SquareMatrix3 jacCylToCart = Acts::SquareMatrix3::Zero();
  jacCylToCart(0, 0) = std::cos(meanPhi);
  jacCylToCart(0, 1) = -meanRho * std::sin(meanPhi);

  jacCylToCart(1, 0) = std::sin(meanPhi);
  jacCylToCart(1, 1) = meanRho * std::cos(meanPhi);

  jacCylToCart(2, 2) = 1;

  Acts::SquareMatrix3 dirCov = jacCylToCart * covCyl * jacCylToCart.transpose();

  m_cov = Acts::SquareMatrix4::Zero();
  m_cov.block(0, 0, 3, 3) = dirCov;
  m_cov(3, 3) = varP;

  double meanX = meanP * meanRho * std::cos(meanPhi);
  double meanY = meanP * meanRho * std::sin(meanPhi);
  meanZ *= meanP;
  m_mean = Acts::Vector3(meanX, meanY, meanZ);
}

Acts::Vector3 SphericalMomentumGenerator::genMomentum(RandomEngine& rng) const {
  std::uniform_real_distribution<double> uniform(0, 1);

  double pMag = m_cfg.pRange.first +
                (m_cfg.pRange.second - m_cfg.pRange.first) * uniform(rng);

  double phi = m_cfg.phiRange.first +
               (m_cfg.phiRange.second - m_cfg.phiRange.first) * uniform(rng);

  double theta = std::acos(
      std::cos(m_cfg.thetaRange.first) -
      (std::cos(m_cfg.thetaRange.first) - std::cos(m_cfg.thetaRange.second)) *
          uniform(rng));

  return pMag * Acts::Vector3(std::sin(theta) * std::cos(phi),
                              std::sin(theta) * std::sin(phi), std::cos(theta));
}

Acts::SquareMatrix4 SphericalMomentumGenerator::getCovariance() const {
  return m_cov;
}

Acts::Vector3 SphericalMomentumGenerator::getMean() const {
  return m_mean;
}
