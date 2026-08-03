#include "TrackingPipeline/TrackFinding/E320TrackParametersEstimator.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "TrackingPipeline/Geometry/E320GeometryOptions.hpp"
#include "TrackingPipeline/MagneticField/ConstantMagField.hpp"
#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"
#include "TrackingPipeline/TrackFitting/StraightLineGX2Fitter.hpp"
#include "TrackingPipeline/Utilities/ThetaMcsRmsCalculator.hpp"

using namespace Acts::UnitLiterals;

namespace E320 {

E320TrackParametersEstimator::E320TrackParametersEstimator(const Config& config)
    : m_cfg(config) {
  const auto& goInst = *E320::GeometryOptions::instance();

  m_dipoleDirIdx = goInst.dipoleDirIdx;
  m_dipoleId = goInst.dipoleId;
  m_dipoleLength = goInst.dipoleHalfPrimary * 2;
  m_dipoleFieldStrength = goInst.dipoleFieldStrength;
  m_sensorThickness = goInst.pixelThickness;
  m_particleCharge =
      (m_cfg.propDirection == PropagationDirection::forward) ? 1_e : -1_e;
}

Acts::BoundMatrix E320TrackParametersEstimator::transportCovToReference(
    const Acts::GeometryContext& gctx, const Acts::Vector3& refSurfacePoint,
    const Acts::Vector3& point, const Acts::Vector3& dir,
    const StraightLineGX2Fitter::Covariance& cov,
    double dipoleFieldStrength) const {
  using enum StraightLineGX2Fitter::GX2FreeIndices;

  const Acts::Vector3& refSurfCenter = m_cfg.referenceSurface->center(gctx);
  const Acts::Vector3& refSurfNormal = m_cfg.referenceSurface->normal(
      gctx, refSurfCenter, Acts::Vector3::UnitX());

  double dirDotProd = dir.dot(refSurfNormal);
  double posDotProd = (refSurfCenter - point).dot(refSurfNormal);

  // Transport covariance to the reference surface
  // and add qOverP and time dimensions
  Acts::ActsMatrix<Acts::eFreeSize, eFreeSize> transportJac;
  transportJac.setZero();

  // Row indices are defined by the Acts
  // Column indices are defined by the GX2 fitter
  Acts::ActsVector<eFreeSize> yPrimeRow;
  yPrimeRow.setZero();
  yPrimeRow(eFreePos0) = -refSurfNormal.x() * dir.y() / dirDotProd;
  yPrimeRow(eFreePos1) = 1 - refSurfNormal.y() * dir.y() / dirDotProd;
  yPrimeRow(eFreePos2) = -refSurfNormal.z() * dir.y() / dirDotProd;
  yPrimeRow(eFreeDir0) =
      -dir.y() * refSurfNormal.x() * posDotProd / (dirDotProd * dirDotProd);
  yPrimeRow(eFreeDir1) =
      posDotProd / dirDotProd -
      dir.y() * refSurfNormal.y() * posDotProd / (dirDotProd * dirDotProd);
  yPrimeRow(eFreeDir2) =
      -dir.y() * refSurfNormal.z() * posDotProd / (dirDotProd * dirDotProd);
  transportJac.block(Acts::eFreePos1, eFreePos0, 1, eFreeSize) =
      yPrimeRow.transpose();

  Acts::ActsVector<eFreeSize> zPrimeRow;
  zPrimeRow.setZero();
  zPrimeRow(eFreePos0) = -refSurfNormal.x() * dir.z() / dirDotProd;
  zPrimeRow(eFreePos1) = -refSurfNormal.y() * dir.z() / dirDotProd;
  zPrimeRow(eFreePos2) = 1 - refSurfNormal.z() * dir.z() / dirDotProd;
  zPrimeRow(eFreeDir0) =
      -dir.z() * refSurfNormal.x() * posDotProd / (dirDotProd * dirDotProd);
  zPrimeRow(eFreeDir1) =
      -dir.z() * refSurfNormal.y() * posDotProd / (dirDotProd * dirDotProd);
  zPrimeRow(eFreeDir2) =
      posDotProd / dirDotProd -
      dir.z() * refSurfNormal.z() * posDotProd / (dirDotProd * dirDotProd);
  transportJac.block(Acts::eFreePos2, eFreePos0, 1, eFreeSize) =
      zPrimeRow.transpose();

  transportJac(Acts::eFreeDir0, eFreeDir0) = 1;
  transportJac(Acts::eFreeDir1, eFreeDir1) = 1;
  transportJac(Acts::eFreeDir2, eFreeDir2) = 1;

  Acts::ActsVector<eFreeSize> pInvRow;
  double pInvDenom = dipoleFieldStrength * m_dipoleLength *
                     std::hypot(dir.x(), dir.y()) / m_particleCharge;
  pInvRow.setZero();
  pInvRow(eFreeDir0) = -dir.y() / pInvDenom;
  pInvRow(eFreeDir1) = dir.x() / pInvDenom;
  transportJac.block(Acts::eFreeQOverP, eFreePos0, 1, eFreeSize) =
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
    const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
    const std::vector<Acts::SourceLink>& sourceLinks,
    const std::vector<std::size_t>& sourceLinkIndices, const Acts::Vector3& dir,
    const Acts::Vector3& point) const {
  Acts::Vector3 newDir = dir;
  Acts::Vector3 newPoint = point;
  StraightLineGX2Fitter::Covariance newCov;

  double chi2 = 0;
  double thetaY = std::atan(newDir.y() / newDir.x());
  double dipoleFieldStrength = m_dipoleFieldStrength;
  if (mctx.hasValue()) {
    const auto& store = mctx.get<std::shared_ptr<MagneticFieldStore>&>();
    if (store->store.contains(m_dipoleId)) {
      dipoleFieldStrength = store->store.at(m_dipoleId)
                                .as<ConstantMagField::Cache>()
                                .m_field(m_dipoleDirIdx);
    }
  }
  double absMom = std::abs(dipoleFieldStrength * m_dipoleLength / thetaY);
  for (std::size_t l = 0; l < m_cfg.nIterations; l++) {
    double thetaRms = detail::getMcpThetaRmsSi(m_sensorThickness, absMom,
                                               std::abs(m_particleCharge));

    chi2 = m_cfg.gx2Fitter->gx2Fit(gctx, sourceLinks, sourceLinkIndices,
                                   thetaRms, newPoint, newDir, newCov);

    thetaY = std::atan(newDir.y() / newDir.x());
    absMom = std::abs(dipoleFieldStrength * m_dipoleLength / thetaY);
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

  dipoleFieldStrength = m_dipoleFieldStrength;
  if (mctx.hasValue()) {
    const auto& store = mctx.get<std::shared_ptr<MagneticFieldStore>&>();
    if (store->store.contains(m_dipoleId)) {
      dipoleFieldStrength = store->store.at(m_dipoleId)
                                .as<ConstantMagField::Cache>()
                                .m_field(m_dipoleDirIdx);
    }
  }
  absMom = std::abs(dipoleFieldStrength * m_dipoleLength / thetaY);

  // Estimate the abs momentum variance
  Acts::BoundMatrix cov = m_cfg.originCov;
  cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = transportCovToReference(
      gctx, vertex3, newPoint, newDir, newCov, dipoleFieldStrength)(
      Acts::eBoundQOverP, Acts::eBoundQOverP);
  if (m_cfg.propDirection == PropagationDirection::forward) {
    return Acts::CurvilinearTrackParameters(
        vertex, newDir, m_particleCharge / absMom, cov,
        Acts::ParticleHypothesis::electron());
  } else {
    return Acts::CurvilinearTrackParameters(
        vertex, -newDir, m_particleCharge / absMom, cov,
        Acts::ParticleHypothesis::electron());
  }
}

}  // namespace E320
