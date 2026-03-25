#include "TrackingPipeline/TrackFinding/E320SeedingAlgorithm.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>

#include <cstddef>
#include <functional>
#include <span>
#include <utility>
#include <vector>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"

using namespace Acts::UnitLiterals;

E320SeedingAlgorithm::E320SeedingAlgorithm(const Config& config,
                                           Acts::Logging::Level level)
    : IAlgorithm("E320SeedingAlgorithm", level), m_cfg(config) {
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_outputSeeds.initialize(m_cfg.outputSeeds);

  const auto& goInst = *E320::GeometryOptions::instance();

  m_dipoleLength = goInst.dipoleHalfPrimary * 2;
  m_dipoleFieldStrength = goInst.dipoleFieldStrength;
}

Acts::BoundMatrix E320SeedingAlgorithm::transportCovToReference(
    const Acts::GeometryContext& gctx, const Acts::Vector3& refSurfacePoint,
    const Acts::Vector3& point, const Acts::Vector3& dir,
    const Acts::ActsSquareMatrix<6>& cov) const {
  const Acts::Vector3& refSurfCenter = m_cfg.referenceSurface->center(gctx);

  // Transport covariance to the reference surface
  // and add qOverP and time dimensions
  Acts::ActsMatrix<8, 6> transportJac;
  transportJac.setZero();

  Acts::ActsVector<6> yPrimeRow;
  yPrimeRow.setZero();
  yPrimeRow(Acts::eFreePos0) = -dir.y() / dir.x();
  yPrimeRow(Acts::eFreePos1) = 1;
  yPrimeRow(Acts::eFreeDir0) =
      -dir.y() * (refSurfCenter.x() - point.x()) / (dir.x() * dir.x());
  yPrimeRow(Acts::eFreeDir1) = (refSurfCenter.x() - point.x()) / dir.x();

  Acts::ActsVector<6> zPrimeRow;
  zPrimeRow.setZero();
  zPrimeRow(Acts::eFreePos0) = -dir.z() / dir.x();
  zPrimeRow(Acts::eFreePos2) = 1;
  zPrimeRow(Acts::eFreeDir0) =
      -dir.z() * (refSurfCenter.x() - point.x()) / (dir.x() * dir.x());
  zPrimeRow(Acts::eFreeDir2) = (refSurfCenter.x() - point.x()) / dir.x();

  Acts::ActsVector<6> pInvRow;
  double pInvDenom = m_dipoleFieldStrength * m_dipoleLength *
                     std::pow(dir.x() * dir.x() + dir.y() * dir.y(), 1.5);
  pInvRow.setZero();
  pInvRow(Acts::eFreeDir0) = -dir.y() * dir.x() / pInvDenom;
  pInvRow(Acts::eFreeDir1) = dir.x() * dir.x() / pInvDenom;

  transportJac.block(Acts::eFreePos1, Acts::eFreePos0, 1, 6) =
      yPrimeRow.transpose();
  transportJac.block(Acts::eFreePos2, Acts::eFreePos0, 1, 6) =
      zPrimeRow.transpose();

  transportJac(Acts::eFreeDir0, Acts::eFreeDir0) = 1;
  transportJac(Acts::eFreeDir1, Acts::eFreeDir1) = 1;
  transportJac(Acts::eFreeDir2, Acts::eFreeDir2) = 1;

  transportJac.block(Acts::eFreeQOverP, Acts::eFreePos0, 1, 6) =
      pInvRow.transpose();

  Acts::FreeMatrix transportedCov =
      transportJac * cov * transportJac.transpose();
  transportedCov(Acts::eFreeTime, Acts::eFreeTime) = 25_ns;

  Acts::FreeToBoundMatrix jacToLoc =
      m_cfg.referenceSurface->freeToBoundJacobian(gctx, refSurfacePoint, dir);
  return jacToLoc * transportedCov * jacToLoc.transpose();
};

ProcessCode E320SeedingAlgorithm::execute(const AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

  const auto& inputSourceLinks = m_inputSourceLinks(ctx);

  ACTS_DEBUG("Received " << inputSourceLinks.size() << " source links");

  if (inputSourceLinks.empty()) {
    ACTS_DEBUG("Input is empty. Skipping");
    m_outputSeeds(ctx, Seeds());
    return ProcessCode::SUCCESS;
  }

  Seeds outSeeds;
  std::vector<Acts::SourceLink> bpmSourceLinks;
  std::vector<std::reference_wrapper<const Acts::SourceLink>> sourceLinkRefs;
  sourceLinkRefs.reserve(inputSourceLinks.size());
  for (const auto& sl : inputSourceLinks) {
    int geoId = sl.get<SimpleSourceLink>().geometryId().sensitive();
    if (geoId >= m_cfg.htOptions.firstLayerId &&
        geoId <= m_cfg.htOptions.lastLayerId) {
      sourceLinkRefs.push_back(std::cref(sl));
    }
    if (std::find(m_cfg.bpmIds.begin(), m_cfg.bpmIds.end(), geoId) !=
        m_cfg.bpmIds.end()) {
      bpmSourceLinks.push_back(sl);
    }
  }
  sourceLinkRefs.shrink_to_fit();

  std::vector<HoughTransformSeeder::HTSeed> htSeeds = m_cfg.htSeeder->findSeeds(
      ctx.geoContext, sourceLinkRefs, m_cfg.htOptions);

  const Acts::Vector3& refSurfCenter =
      m_cfg.referenceSurface->center(ctx.geoContext);
  const Acts::Vector3& refSurfNormal = m_cfg.referenceSurface->normal(
      ctx.geoContext, refSurfCenter, Acts::Vector3::UnitX());

  ACTS_DEBUG("Found " << htSeeds.size() << " HT seeds");
  outSeeds.reserve(htSeeds.size());
  for (std::size_t i = 0; i < htSeeds.size(); i++) {
    const auto& [point, dir, cov, sl, chi2, count] = htSeeds.at(i);

    // Transport parameteres to the reference surface
    double dVertex =
        (refSurfCenter - point).dot(refSurfNormal) / dir.dot(refSurfNormal);
    Acts::Vector3 vertex3 = point + dir * dVertex;
    Acts::Vector4 vertex(vertex3.x(), vertex3.y(), vertex3.z(), 0);

    double thetaY = std::atan(dir.y() / dir.x());
    double absMom =
        std::abs(m_dipoleFieldStrength * m_dipoleLength / std::sin(thetaY));

    // Estimate the abs momentum variance
    Acts::BoundMatrix originCov = m_cfg.originCov;
    // originCov(Acts::eBoundQOverP, Acts::eBoundQOverP) =
    //     transportCovToReference(ctx.geoContext, vertex3, point, dir, cov)(
    //         Acts::eBoundQOverP, Acts::eBoundQOverP);

    std::vector<Acts::SourceLink> sourceLinks = sl;
    sourceLinks.insert(sourceLinks.end(), bpmSourceLinks.begin(),
                       bpmSourceLinks.end());

    if (m_cfg.propDirection == PropagationDirection::forward) {
      outSeeds.emplace_back(sourceLinks,
                            Acts::CurvilinearTrackParameters(
                                vertex, dir, -1_e / absMom, originCov,
                                Acts::ParticleHypothesis::electron()),
                            i);
    } else if (m_cfg.propDirection == PropagationDirection::backward) {
      outSeeds.emplace_back(sourceLinks,
                            Acts::CurvilinearTrackParameters(
                                vertex, -dir, 1_e / absMom, originCov,
                                Acts::ParticleHypothesis::electron()),
                            i);
    }
  }
  ACTS_DEBUG("Found " << outSeeds.size() << " seeds");
  ACTS_DEBUG("Sending " << outSeeds.size() << " seeds");
  m_outputSeeds(ctx, std::move(outSeeds));
  return ProcessCode::SUCCESS;
}
