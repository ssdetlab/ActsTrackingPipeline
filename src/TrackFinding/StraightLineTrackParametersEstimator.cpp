#include "TrackingPipeline/TrackFinding/StraightLineTrackParametersEstimator.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cstddef>
#include <vector>

StraightLineTrackParametersEstimator::StraightLineTrackParametersEstimator(
    const Config& config)
    : m_cfg(config) {}

StraightLineTrackParametersEstimator::Result
StraightLineTrackParametersEstimator::estimateParameters(
    const Acts::GeometryContext& gctx,
    const std::vector<Acts::SourceLink>& sourceLinks,
    const std::vector<std::size_t>& sourceLinkIndices, const Acts::Vector3& dir,
    const Acts::Vector3& point) const {
  Acts::Vector3 newDir = dir;
  Acts::Vector3 newPoint = point;
  StraightLineGX2Fitter::Covariance newCov;

  double chi2 = 0;
  for (std::size_t l = 0; l < m_cfg.nIterations; l++) {
    chi2 = m_cfg.gx2Fitter->gx2Fit(gctx, sourceLinks, sourceLinkIndices,
                                   m_cfg.thetaRms, newPoint, newDir, newCov);
  }
  if (chi2 > m_cfg.maxChi2) {
    return Acts::Result<Acts::CurvilinearTrackParameters>::failure(
        std::error_code());
  }

  const Acts::Vector3& refSurfCenter = m_cfg.referenceSurface->center(gctx);
  const Acts::Vector3& refSurfNormal = m_cfg.referenceSurface->normal(
      gctx, refSurfCenter, Acts::Vector3::UnitX());

  std::cout << "POINT " << newPoint.transpose() << "\n";
  std::cout << "DIR " << newDir.transpose() << "\n";

  // Transport parameteres to the reference surface
  double dVertex =
      (refSurfCenter - newPoint).dot(refSurfNormal) / newDir.dot(refSurfNormal);
  Acts::Vector3 vertex3 = newPoint + newDir * dVertex;
  Acts::Vector4 vertex(vertex3.x(), vertex3.y(), vertex3.z(), 0);
  std::cout << "VERTEX " << vertex.transpose() << "\n";

  return Acts::CurvilinearTrackParameters(
      vertex, newDir, m_cfg.charge / m_cfg.absMomentum, m_cfg.originCov,
      m_cfg.particleHypothesis);
}
