#include "TrackingPipeline/Alignment/LocalAlignmentParametersSolverSVD.hpp"

#include "Acts/Utilities/Logger.hpp"

LocalAlignmentParametersSolverSVD::LocalAlignmentParametersSolverSVD(
    const Config& cfg, Acts::Logging::Level level)
    : m_cfg(cfg),
      m_logger(
          Acts::getDefaultLogger("LocalAlignmentParametersSolverSVD", level)) {
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
};

void LocalAlignmentParametersSolverSVD::calculateAlignmentParameters(
    const Acts::GeometryContext& /*gctx*/,
    ActsAlignment::AlignmentResult& alignRes,
    const Acts::ActsDynamicVector& sumChi2Derivative,
    const Acts::ActsDynamicMatrix& sumChi2SecondDerivative) const {
  Eigen::JacobiSVD<Acts::ActsDynamicMatrix> svd;
  svd.compute(sumChi2SecondDerivative,
              Eigen::ComputeThinU | Eigen::ComputeThinV);

  Acts::ActsDynamicVector singularValues = svd.singularValues();
  ACTS_INFO("Hessian singular values:\n" << singularValues);

  double maxSv = singularValues.maxCoeff();
  Acts::ActsDynamicVector invSingularValues(singularValues.size());

  std::size_t cutoffIdx = 1;
  for (std::size_t i = 1; i < singularValues.size(); i++) {
    if (singularValues(i) < maxSv * m_cfg.maxSingularValueTol &&
        singularValues(i) / singularValues(i - 1) < m_cfg.singularValueGapTol) {
      cutoffIdx = i;
      break;
    }
  }

  invSingularValues(0) = 1.0 / singularValues(0);
  for (std::size_t i = 1; i < singularValues.size(); i++) {
    double sv = singularValues(i);
    if (i < cutoffIdx) {
      invSingularValues(i) = 1.0 / sv;
      ACTS_INFO("Accepting mode with singular value "
                << sv << " and ratios: " << sv / maxSv << ", "
                << sv / singularValues(i - 1));
    } else {
      invSingularValues(i) = 0.0;
      ACTS_INFO("Suppressing weak mode " << i << " with singular value " << sv
                                         << " and ratios: " << sv / maxSv
                                         << ", " << sv / singularValues(i - 1));
    }
  }

  Acts::ActsDynamicMatrix pinv = svd.matrixV() *
                                 invSingularValues.asDiagonal() *
                                 svd.matrixU().transpose();

  alignRes.deltaAlignmentParameters = -pinv * sumChi2Derivative;

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
  // Alignment parameters covariance
  alignRes.alignmentCovariance = 2 * sumChi2SecondDerivativeActive.inverse();
  alignRes.alignmentDof = activeDofs;

  // Chi2 change
  alignRes.deltaChi2 =
      0.5 * sumChi2Derivative.dot(alignRes.deltaAlignmentParameters);

  ACTS_INFO("Alignment parameters covariance:\n"
            << alignRes.alignmentCovariance);
  ACTS_INFO("Number of degrees of freedom " << alignRes.alignmentDof);
  ACTS_INFO("Delta chi2 " << alignRes.deltaChi2);
}
