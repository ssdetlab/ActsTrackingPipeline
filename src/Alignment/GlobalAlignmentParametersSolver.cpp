#include "TrackingPipeline/Alignment/GlobalAlignmentParametersSolver.hpp"

GlobalAlignmentParametersSolver::GlobalAlignmentParametersSolver(
    const Config& cfg, Acts::Logging::Level level)
    : m_cfg(cfg),
      m_logger(
          Acts::getDefaultLogger("GlobalAlignmentParametersSolver", level)) {
  m_doCenter0 = ACTS_CHECK_BIT(m_cfg.alignmentMask,
                               ActsAlignment::AlignmentMask::Center0);
  m_doCenter1 = ACTS_CHECK_BIT(m_cfg.alignmentMask,
                               ActsAlignment::AlignmentMask::Center1);
  m_doCenter2 = ACTS_CHECK_BIT(m_cfg.alignmentMask,
                               ActsAlignment::AlignmentMask::Center2);
  m_doAngle0 = ACTS_CHECK_BIT(m_cfg.alignmentMask,
                              ActsAlignment::AlignmentMask::Rotation0);
  m_doAngle1 = ACTS_CHECK_BIT(m_cfg.alignmentMask,
                              ActsAlignment::AlignmentMask::Rotation1);
  m_doAngle2 = ACTS_CHECK_BIT(m_cfg.alignmentMask,
                              ActsAlignment::AlignmentMask::Rotation2);

  if (m_doCenter0) {
    m_activeIdxs.push_back(Acts::eAlignmentCenter0);
  }
  if (m_doCenter1) {
    m_activeIdxs.push_back(Acts::eAlignmentCenter1);
  }
  if (m_doCenter2) {
    m_activeIdxs.push_back(Acts::eAlignmentCenter2);
  }
  if (m_doAngle0) {
    m_activeIdxs.push_back(Acts::eAlignmentRotation0);
  }
  if (m_doAngle1) {
    m_activeIdxs.push_back(Acts::eAlignmentRotation1);
  }
  if (m_doAngle2) {
    m_activeIdxs.push_back(Acts::eAlignmentRotation2);
  }
};

void GlobalAlignmentParametersSolver::calculateAlignmentParameters(
    const Acts::GeometryContext& gctx, ActsAlignment::AlignmentResult& alignRes,
    const Acts::ActsDynamicVector& sumChi2Derivative,
    const Acts::ActsDynamicMatrix& sumChi2SecondDerivative) const {
  std::size_t nAlignDof = alignRes.alignmentDof;

  Acts::ActsDynamicMatrix U =
      Acts::ActsDynamicMatrix::Zero(nAlignDof, m_nRigidDof);

  Acts::Vector3 centerOfMass = Acts::Vector3::Zero();
  for (const auto& [surf, idx] : alignRes.idxedAlignSurfaces) {
    centerOfMass += surf->center(gctx);
  }
  centerOfMass /= alignRes.idxedAlignSurfaces.size();

  for (const auto& [surf, idx] : alignRes.idxedAlignSurfaces) {
    if (m_doAngle0) {
      U.block(idx * Acts::eAlignmentSize + Acts::eAlignmentCenter0, 3, 3, 1) =
          Acts::Vector3::UnitX().cross(surf->center(gctx) - centerOfMass);
      U(idx * Acts::eAlignmentSize + Acts::eAlignmentRotation0, 3) = 1;
    }
    if (m_doAngle1) {
      U.block(idx * Acts::eAlignmentSize + Acts::eAlignmentCenter0, 4, 3, 1) =
          Acts::Vector3::UnitY().cross(surf->center(gctx) - centerOfMass);
      U(idx * Acts::eAlignmentSize + Acts::eAlignmentRotation1, 4) = 1;
    }
    if (m_doAngle2) {
      U.block(idx * Acts::eAlignmentSize + Acts::eAlignmentCenter0, 5, 3, 1) =
          Acts::Vector3::UnitZ().cross(surf->center(gctx) - centerOfMass);
      U(idx * Acts::eAlignmentSize + Acts::eAlignmentRotation2, 5) = 1;
    }

    if (m_doCenter0) {
      U(idx * Acts::eAlignmentSize + Acts::eAlignmentCenter0, 0) = 1;
    }
    if (m_doCenter1) {
      U(idx * Acts::eAlignmentSize + Acts::eAlignmentCenter1, 1) = 1;
    }
    if (m_doCenter2) {
      U(idx * Acts::eAlignmentSize + Acts::eAlignmentCenter2, 2) = 1;
    }
  }

  Acts::ActsDynamicMatrix rigidA = U.transpose() * sumChi2SecondDerivative * U;
  Acts::ActsDynamicVector rigidb = -U.transpose() * sumChi2Derivative;

  double eps = 1e-12 * std::max(1.0, rigidA.norm());
  if (rigidA.fullPivLu().rank() < m_nRigidDof) {
    rigidA.diagonal().array() += eps;
  }

  Acts::ActsDynamicVector rigidDelta = rigidA.householderQr().solve(rigidb);
  alignRes.deltaAlignmentParameters = U * rigidDelta;

  ACTS_INFO("Delta alignment parameters\n"
            << alignRes.deltaAlignmentParameters);

  // Alignment parameters covariance
  std::size_t activeDofs = m_activeIdxs.size();
  Acts::ActsDynamicMatrix rigidAActive =
      Acts::ActsDynamicMatrix::Zero(activeDofs, activeDofs);
  for (std::size_t i = 0; i < m_activeIdxs.size(); i++) {
    std::size_t idxA = m_activeIdxs.at(i);
    for (std::size_t j = 0; j < m_activeIdxs.size(); j++) {
      std::size_t idxB = m_activeIdxs.at(j);
      rigidAActive(i, j) = rigidA(idxA, idxB);
    }
  }

  alignRes.alignmentCovariance = 2 * rigidAActive.inverse();
  alignRes.alignmentDof = activeDofs;

  // Chi2 change
  alignRes.deltaChi2 =
      0.5 * sumChi2Derivative.dot(alignRes.deltaAlignmentParameters);
}
