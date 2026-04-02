#include "TrackingPipeline/Alignment/LocalAlignmentParametersSolverConstraints.hpp"

#include <array>

LocalAlignmentParametersSolverConstraints::
    LocalAlignmentParametersSolverConstraints(const Config& cfg,
                                              Acts::Logging::Level level)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger(
          "LocalAlignmentParametersSolverConstraints", level)) {
  std::array<ActsAlignment::AlignmentMask, 6> maskElements{
      ActsAlignment::AlignmentMask::Center0,
      ActsAlignment::AlignmentMask::Center1,
      ActsAlignment::AlignmentMask::Center2,
      ActsAlignment::AlignmentMask::Rotation0,
      ActsAlignment::AlignmentMask::Rotation1,
      ActsAlignment::AlignmentMask::Rotation2};
  std::array<std::size_t, 6> alignmentDofs{
      Acts::eAlignmentCenter0,   Acts::eAlignmentCenter1,
      Acts::eAlignmentCenter2,   Acts::eAlignmentRotation0,
      Acts::eAlignmentRotation1, Acts::eAlignmentRotation2};

  for (std::size_t i = 0; i < maskElements.size(); i++) {
    if (ACTS_CHECK_BIT(m_cfg.alignmentMask, maskElements.at(i))) {
      m_activeIdxs.push_back(alignmentDofs.at(i));
    }
  }
}

void LocalAlignmentParametersSolverConstraints::calculateAlignmentParameters(
    const Acts::GeometryContext& /*gctx*/,
    ActsAlignment::AlignmentResult& alignRes,
    const Acts::ActsDynamicVector& sumChi2Derivative,
    const Acts::ActsDynamicMatrix& sumChi2SecondDerivative) const {
  std::size_t nAlignDof = alignRes.alignmentDof;

  ACTS_INFO("Unconstrained system rank = "
            << sumChi2SecondDerivative.fullPivHouseholderQr().rank() << "/"
            << sumChi2SecondDerivative.rows());

  // Constraint handling
  Acts::ActsDynamicMatrix cMat =
      Acts::ActsDynamicMatrix::Zero(m_nConstraints, nAlignDof);
  Acts::ActsDynamicVector cVec = Acts::ActsDynamicVector::Zero(m_nConstraints);

  for (const auto& [surf, idx] : alignRes.idxedAlignSurfaces) {
    cMat(0, idx * Acts::eAlignmentSize + Acts::eAlignmentCenter1) = 1.0;
    cMat(1, idx * Acts::eAlignmentSize + Acts::eAlignmentCenter2) = 1.0;
    cMat(2, idx * Acts::eAlignmentSize + Acts::eAlignmentRotation2) = 1.0;
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
  ACTS_INFO("Delta alignment parameters sum\n"
            << cMat * alignRes.deltaAlignmentParameters);

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
