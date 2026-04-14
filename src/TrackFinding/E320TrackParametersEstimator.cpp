#include "TrackingPipeline/TrackFinding/E320TrackParametersEstimator.hpp"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/SourceLink.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Utilities/Result.hpp>

#include <cstddef>
#include <vector>

#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Utilities/ThetaMcsRmsCalculator.hpp"

using namespace Acts::UnitLiterals;

namespace E320 {

E320TrackParametersEstimator::E320TrackParametersEstimator(const Config& config)
    : m_cfg(config) {
  const auto& goInst = *E320::GeometryOptions::instance();

  m_dipoleLength = goInst.dipoleHalfPrimary * 2;
  m_dipoleFieldStrength = goInst.dipoleFieldStrength;
}

Acts::BoundMatrix E320TrackParametersEstimator::transportCovToReference(
    const Acts::GeometryContext& gctx, const Acts::Vector3& refSurfacePoint,
    const Acts::Vector3& point, const Acts::Vector3& dir,
    const StraightLineGX2Fitter::Covariance& cov) const {
  constexpr std::size_t gx2CovDim = StraightLineGX2Fitter::covarianceDim;

  const Acts::Vector3& refSurfCenter = m_cfg.referenceSurface->center(gctx);

  // Transport covariance to the reference surface
  // and add qOverP and time dimensions
  Acts::ActsMatrix<8, gx2CovDim> transportJac;
  transportJac.setZero();

  Acts::ActsVector<gx2CovDim> yPrimeRow;
  yPrimeRow.setZero();
  yPrimeRow(Acts::eFreePos0) = -dir.y() / dir.x();
  yPrimeRow(Acts::eFreePos1) = 1;
  yPrimeRow(Acts::eFreeDir0) =
      -dir.y() * (refSurfCenter.x() - point.x()) / (dir.x() * dir.x());
  yPrimeRow(Acts::eFreeDir1) = (refSurfCenter.x() - point.x()) / dir.x();

  Acts::ActsVector<gx2CovDim> zPrimeRow;
  zPrimeRow.setZero();
  zPrimeRow(Acts::eFreePos0) = -dir.z() / dir.x();
  zPrimeRow(Acts::eFreePos2) = 1;
  zPrimeRow(Acts::eFreeDir0) =
      -dir.z() * (refSurfCenter.x() - point.x()) / (dir.x() * dir.x());
  zPrimeRow(Acts::eFreeDir2) = (refSurfCenter.x() - point.x()) / dir.x();

  Acts::ActsVector<gx2CovDim> pInvRow;
  double pInvDenom = m_dipoleFieldStrength * m_dipoleLength *
                     std::pow(dir.x() * dir.x() + dir.y() * dir.y(), 1.5);
  pInvRow.setZero();
  pInvRow(Acts::eFreeDir0) = -dir.y() * dir.x() / pInvDenom;
  pInvRow(Acts::eFreeDir1) = dir.x() * dir.x() / pInvDenom;

  transportJac.block(Acts::eFreePos1, Acts::eFreePos0, 1, gx2CovDim) =
      yPrimeRow.transpose();
  transportJac.block(Acts::eFreePos2, Acts::eFreePos0, 1, gx2CovDim) =
      zPrimeRow.transpose();

  transportJac(Acts::eFreeDir0, Acts::eFreeDir0) = 1;
  transportJac(Acts::eFreeDir1, Acts::eFreeDir1) = 1;
  transportJac(Acts::eFreeDir2, Acts::eFreeDir2) = 1;

  transportJac.block(Acts::eFreeQOverP, Acts::eFreePos0, 1, gx2CovDim) =
      pInvRow.transpose();

  Acts::FreeMatrix transportedCov =
      transportJac * cov * transportJac.transpose();
  transportedCov(Acts::eFreeTime, Acts::eFreeTime) = 1_fs;

  Acts::FreeToBoundMatrix jacToLoc =
      m_cfg.referenceSurface->freeToBoundJacobian(gctx, refSurfacePoint, dir);
  return jacToLoc * transportedCov * jacToLoc.transpose();
};

E320TrackParametersEstimator::Result
E320TrackParametersEstimator::estimateParameters(
    const Acts::GeometryContext& gctx,
    const std::vector<Acts::SourceLink>& sourceLinks,
    const std::vector<std::size_t>& sourceLinkIndices, const Acts::Vector3& dir,
    const Acts::Vector3& point) const {
  Acts::Vector3 newDir = dir;
  Acts::Vector3 newPoint = point;
  StraightLineGX2Fitter::Covariance newCov;

  double chi2 = 0;
  double thetaY = std::atan(newDir.y() / newDir.x());
  double absMom = std::abs(m_dipoleFieldStrength * m_dipoleLength / thetaY);
  for (std::size_t l = 0; l < m_cfg.nIterations; l++) {
    double sensorThickness = 50_um;
    double particleCharge = 1_e;

    double thetaRms =
        detail::getMcpThetaRmsSi(sensorThickness, absMom, particleCharge);

    chi2 = m_cfg.gx2Fitter->gx2Fit(gctx, sourceLinks, sourceLinkIndices,
                                   thetaRms, newPoint, newDir, newCov);
  }
  if (chi2 > m_cfg.maxChi2) {
    return Acts::Result<Acts::CurvilinearTrackParameters>::failure(
        std::error_code());
  }

  const Acts::Vector3& refSurfCenter = m_cfg.referenceSurface->center(gctx);
  const Acts::Vector3& refSurfNormal = m_cfg.referenceSurface->normal(
      gctx, refSurfCenter, Acts::Vector3::UnitX());

  // Transport parameteres to the reference surface
  double dVertex =
      (refSurfCenter - newPoint).dot(refSurfNormal) / newDir.dot(refSurfNormal);
  Acts::Vector3 vertex3 = newPoint + newDir * dVertex;
  Acts::Vector4 vertex(vertex3.x(), vertex3.y(), vertex3.z(), 0);

  thetaY = std::atan(newDir.y() / newDir.x());
  absMom = std::abs(m_dipoleFieldStrength * m_dipoleLength / thetaY);

  // Estimate the abs momentum variance
  Acts::BoundMatrix cov = m_cfg.originCov;
  cov(Acts::eBoundQOverP, Acts::eBoundQOverP) =
      transportCovToReference(gctx, vertex3, newPoint, newDir, newCov)(
          Acts::eBoundQOverP, Acts::eBoundQOverP);

  if (m_cfg.propDirection == PropagationDirection::forward) {
    return Acts::CurvilinearTrackParameters(
        vertex, newDir, 1_e / absMom, cov,
        Acts::ParticleHypothesis::electron());
  } else {
    return Acts::CurvilinearTrackParameters(
        vertex, -newDir, -1_e / absMom, cov,
        Acts::ParticleHypothesis::electron());
  }
}

}  // namespace E320
