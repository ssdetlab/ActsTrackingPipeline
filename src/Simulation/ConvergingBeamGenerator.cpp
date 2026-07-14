#include "TrackingPipeline/Simulation/ConvergingBeamGenerator.hpp"

#include <Acts/Definitions/Algebra.hpp>

#include <cmath>

ConvergingBeamGenerator::ConvergingBeamGenerator(const Config& cfg)
    : m_cfg(cfg), m_state(std::make_unique<State>()) {
  m_state->genState = {false, false};
}

Acts::Vector3 ConvergingBeamGenerator::genVertex(RandomEngine& rng) const {
  if (!m_state->genState.first) {
    internalUpdate(rng);
  }
  m_state->genState.first = false;
  return m_state->vertex;
};

Acts::Vector3 ConvergingBeamGenerator::genMomentum(RandomEngine& rng) const {
  if (!m_state->genState.second) {
    internalUpdate(rng);
  }
  m_state->genState.second = false;
  return m_state->momentum;
};

Acts::SquareMatrix3 ConvergingBeamGenerator::getVertexCovariance() const {
  return Acts::SquareMatrix3::Identity();
}

Acts::SquareMatrix4 ConvergingBeamGenerator::getMomentumCovariance() const {
  return Acts::SquareMatrix4::Identity();
}

Acts::Vector3 ConvergingBeamGenerator::getVertexMean() const {
  return Acts::Vector3::Zero();
}

Acts::Vector3 ConvergingBeamGenerator::getMomentumMean() const {
  return Acts::Vector3::Zero();
}

void ConvergingBeamGenerator::internalUpdate(RandomEngine& rng) const {
  // Sample phase space at the waist
  double waistLong = m_cfg.waistSigmaLong * m_unitNormal(rng);
  double waistShort = m_cfg.waistSigmaShort * m_unitNormal(rng);

  double waistThetaLong = m_cfg.waistSigmaThetaLong * m_unitNormal(rng);
  double waistThetaShort = m_cfg.waistSigmaThetaShort * m_unitNormal(rng);

  double deltaPrimary =
      m_cfg.waistPosition(m_cfg.primaryIdx) - m_cfg.referencePositionPrimary;

  // Transport backwards from waist to reference plane
  double refLong = m_cfg.waistPosition(m_cfg.longIdx) + waistLong -
                   deltaPrimary * waistThetaLong;
  double refShort = m_cfg.waistPosition(m_cfg.shortIdx) + waistShort -
                    deltaPrimary * waistThetaShort;
  double refPrimary = m_cfg.referencePositionPrimary;

  Acts::Vector3 position(0, 0, 0);
  position(m_cfg.primaryIdx) = refPrimary;
  position(m_cfg.longIdx) = refLong;
  position(m_cfg.shortIdx) = refShort;

  // Direction corresponding to slopes tx, ty
  Acts::Vector3 direction(0, 0, 0);
  direction(m_cfg.primaryIdx) = 1.0;
  direction(m_cfg.longIdx) = waistThetaLong;
  direction(m_cfg.shortIdx) = waistThetaShort;
  direction.normalize();

  m_state->vertex = position;
  m_state->momentum = m_cfg.beamEnergy * direction;
  m_state->genState = {true, true};
}
