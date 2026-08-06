#include "TrackingPipeline/Alignment/GlobalAlignmentParametersSolverConstraints.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"

#include <algorithm>
#include <cstddef>

GlobalAlignmentParametersSolverConstraints::
    GlobalAlignmentParametersSolverConstraints(const Config& cfg,
                                               Acts::Logging::Level level)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("GlobalAlignmentParametersSolverTest",
                                      level)) {
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

void GlobalAlignmentParametersSolverConstraints::calculateAlignmentParameters(
    const Acts::GeometryContext& gctx, ActsAlignment::AlignmentResult& alignRes,
    const Acts::ActsDynamicVector& sumChi2Derivative,
    const Acts::ActsDynamicMatrix& sumChi2SecondDerivative) const {
  std::size_t nAlignDof = alignRes.alignmentDof;

  Acts::Vector3 centerOfMass = Acts::Vector3::Zero();
  for (const auto& [surf, idx] : alignRes.idxedAlignSurfaces) {
    centerOfMass += surf->center(gctx);
  }
  centerOfMass /= alignRes.idxedAlignSurfaces.size();

  std::vector<Acts::Vector3> signedDistances;
  for (const auto& [surf, idx] : alignRes.idxedAlignSurfaces) {
    signedDistances.push_back(centerOfMass - surf->center(gctx));
  }
  std::sort(signedDistances.begin(), signedDistances.end(),
            [](const auto& a, const auto& b) { return a.x() < b.x(); });

  // Constraint handling
  std::size_t nAngleConstraints = alignRes.idxedAlignSurfaces.size() - 1;
  std::size_t nXConstraints = alignRes.idxedAlignSurfaces.size();
  std::size_t nYConstraints = alignRes.idxedAlignSurfaces.size();
  std::size_t nCenter2Constraints = alignRes.idxedAlignSurfaces.size() - 1;

  std::size_t nConstraints =
      nAngleConstraints + nXConstraints + nYConstraints + nCenter2Constraints;
  Acts::ActsDynamicMatrix cMat = Acts::ActsDynamicMatrix::Zero(
      nConstraints, Acts::eAlignmentSize * alignRes.idxedAlignSurfaces.size());
  Acts::ActsDynamicVector cVec = Acts::ActsDynamicVector::Zero(nConstraints);

  // Angles are all equal
  for (std::size_t i = 0; i < nAngleConstraints; i++) {
    std::size_t row = i;
    std::size_t colPlus = Acts::eAlignmentRotation2;
    std::size_t colMinus =
        (i + 1) * Acts::eAlignmentSize + Acts::eAlignmentRotation2;
    cMat(row, colPlus) = 1.0;
    cMat(row, colMinus) = -1.0;
  }
  for (std::size_t i = 0; i < nXConstraints; i++) {
    std::size_t row = i + nAngleConstraints;
    std::size_t colPlus = i * Acts::eAlignmentSize + Acts::eAlignmentCenter0;
    std::size_t colMinus = i * Acts::eAlignmentSize + Acts::eAlignmentRotation2;
    cMat(row, colPlus) = 1.0;
    cMat(row, colMinus) = -signedDistances.at(i).y();
  }
  for (std::size_t i = 0; i < nYConstraints; i++) {
    if (i == 2) {
      continue;
    }
    std::size_t row = i + nAngleConstraints + nXConstraints;
    std::size_t colPlus = i * Acts::eAlignmentSize + Acts::eAlignmentCenter1;
    std::size_t colMinus = i * Acts::eAlignmentSize + Acts::eAlignmentRotation2;
    cMat(row, colPlus) = 1.0;
    cMat(row, colMinus) = signedDistances.at(i).x();
    cMat(row, 2 * Acts::eAlignmentSize + Acts::eAlignmentCenter1) = -1.0;
  }
  for (std::size_t i = 0; i < nCenter2Constraints; i++) {
    std::size_t row = i + nAngleConstraints + nXConstraints + nYConstraints;
    std::size_t colPlus = Acts::eAlignmentCenter2;
    std::size_t colMinus =
        (i + 1) * Acts::eAlignmentSize + Acts::eAlignmentCenter2;
    cMat(row, colPlus) = 1.0;
    cMat(row, colMinus) = -1.0;
  }

  std::size_t cRows = cMat.rows();
  std::size_t sumChi2SDRows = sumChi2SecondDerivative.rows();
  std::size_t sumChi2SDCols = sumChi2SecondDerivative.cols();

  // Constrained system of equations
  Acts::ActsDynamicMatrix kktMat = Acts::ActsDynamicMatrix::Zero(
      sumChi2SDRows + cRows, sumChi2SDCols + cRows);
  kktMat.topLeftCorner(sumChi2SDRows, sumChi2SDCols) = sumChi2SecondDerivative;
  kktMat.topRightCorner(sumChi2SDRows, cRows) = cMat.transpose();
  kktMat.bottomLeftCorner(cRows, sumChi2SDCols) = cMat;

  Acts::ActsDynamicVector rhs(sumChi2SDRows + cRows);
  rhs.head(sumChi2SDRows) = -sumChi2Derivative;
  rhs.tail(cRows) = cVec;

  Eigen::HouseholderQR<Acts::ActsDynamicMatrix> hqr(kktMat);
  ACTS_INFO("Constrained system rank = "
            << kktMat.fullPivHouseholderQr().rank() << "/"
            << std::max(kktMat.rows(), kktMat.cols()));

  // Solve
  Acts::ActsDynamicVector sol = hqr.solve(rhs);
  alignRes.deltaAlignmentParameters = sol.head(nAlignDof);

  ACTS_INFO("Delta alignment parameters\n"
            << alignRes.deltaAlignmentParameters);

  // Alignment parameters covariance
  std::size_t activeDofs =
      m_activeIdxs.size() * alignRes.idxedAlignSurfaces.size();
  Acts::ActsDynamicMatrix sumChi2SecondDerivativeActive =
      Acts::ActsDynamicMatrix::Zero(activeDofs, activeDofs);
  for (std::size_t i = 0; i < m_activeIdxs.size(); i++) {
    std::size_t idxA = m_activeIdxs.at(i);
    for (std::size_t j = 0; j < m_activeIdxs.size(); j++) {
      std::size_t idxB = m_activeIdxs.at(j);
      for (const auto& [surfA, surfIdxA] : alignRes.idxedAlignSurfaces) {
        for (const auto& [surfB, surfIdxB] : alignRes.idxedAlignSurfaces) {
          sumChi2SecondDerivativeActive(surfIdxA * m_activeIdxs.size() + i,
                                        surfIdxB * m_activeIdxs.size() + j) =
              sumChi2SecondDerivative(surfIdxA * Acts::eAlignmentSize + idxA,
                                      surfIdxB * Acts::eAlignmentSize + idxB);
        }
      }
    }
  }

  alignRes.alignmentCovariance = 2 * sumChi2SecondDerivativeActive.inverse();
  alignRes.alignmentDof = activeDofs;

  // Chi2 change
  alignRes.deltaChi2 =
      0.5 * sumChi2Derivative.dot(alignRes.deltaAlignmentParameters);
}
