#include "TrackingPipeline/TrackFinding/E320SeedingAlgorithm.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>

#include <cstddef>
#include <utility>
#include <vector>

#include <TH1.h>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"

using namespace Acts::UnitLiterals;

namespace {

void constructTracks(const std::shared_ptr<IdxTree::Node>& root,
                     std::vector<std::size_t>& track,
                     std::vector<std::vector<std::size_t>>& tracks) {
  track.push_back(root->m_idx);
  if (root->children.size() == 0) {
    tracks.push_back(track);
    track.pop_back();
    return;
  }
  for (auto& child : root->children) {
    constructTracks(child, track, tracks);
  }
  track.pop_back();
}

}  // namespace

E320SeedingAlgorithm::E320SeedingAlgorithm(const Config& config,
                                           Acts::Logging::Level level)
    : IAlgorithm("E320SeedingAlgorithm", level), m_cfg(config) {
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_inputDetSourceLinkIndices.initialize(m_cfg.inputDetSourceLinkIndices);
  // m_inputBpmSourceLinkIndices.initialize(m_cfg.inputBpmSourceLinkIndices);

  m_outputSeeds.initialize(m_cfg.outputSeeds);
  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);

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
  transportedCov(Acts::eFreeTime, Acts::eFreeTime) = 1_fs;

  Acts::FreeToBoundMatrix jacToLoc =
      m_cfg.referenceSurface->freeToBoundJacobian(gctx, refSurfacePoint, dir);
  return jacToLoc * transportedCov * jacToLoc.transpose();
};

ProcessCode E320SeedingAlgorithm::execute(const AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

  const auto& inputSourceLinks = m_inputSourceLinks(ctx);
  const auto& inputDetSourceLinksIndices = m_inputDetSourceLinkIndices(ctx);
  const auto& inputBpmSourceLinksIndices =
      (ctx.eventStore.exists(m_cfg.inputBpmSourceLinkIndices))
          ? m_inputBpmSourceLinkIndices(ctx)
          : std::vector<std::size_t>();

  ACTS_DEBUG("Received " << inputSourceLinks.size() << " source links");

  if (inputSourceLinks.empty()) {
    ACTS_DEBUG("Input is empty. Skipping");
    m_outputSeeds(ctx, IndexSeeds());
    m_outputTrackParameters(ctx, {});
    return ProcessCode::SUCCESS;
  }

  IndexSeeds outSeeds;
  std::vector<Acts::CurvilinearTrackParameters> outTrackParameters;
  std::vector<HoughTransformSeeder::HTSeed> htSeeds =
      m_cfg.htSeeder->findSeeds(ctx.geoContext, inputSourceLinks,
                                inputDetSourceLinksIndices, m_cfg.htOptions);

  const Acts::Vector3& refSurfCenter =
      m_cfg.referenceSurface->center(ctx.geoContext);
  const Acts::Vector3& refSurfNormal = m_cfg.referenceSurface->normal(
      ctx.geoContext, refSurfCenter, Acts::Vector3::UnitX());

  ACTS_DEBUG("Found " << htSeeds.size() << " HT seeds");
  outSeeds.reserve(htSeeds.size() * m_cfg.maxLayers);
  outTrackParameters.reserve(htSeeds.size() * m_cfg.maxLayers);
  for (std::size_t i = 0; i < htSeeds.size(); i++) {
    const auto& [point, dir, slIdxs] = htSeeds.at(i);

    IdxTree::IdxContainer idxContainer;
    idxContainer.reserve(slIdxs.size());
    std::size_t firstLayerId = std::numeric_limits<std::size_t>::max();
    for (std::size_t idx : slIdxs) {
      std::size_t geoId = inputSourceLinks.at(idx)
                              .get<SimpleSourceLink>()
                              .geometryId()
                              .sensitive();
      idxContainer.push_back({idx, geoId});
      if (firstLayerId > geoId) {
        firstLayerId = geoId;
      }
    }

    std::sort(idxContainer.begin(), idxContainer.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
    auto rootEndIt = std::find_if(
        idxContainer.begin(), idxContainer.end(),
        [&firstLayerId](const auto& a) { return (a.second != firstLayerId); });

    for (auto it = idxContainer.begin(); it != rootEndIt; it++) {
      std::vector<std::size_t> trackContainer;
      trackContainer.reserve(m_cfg.minLayers);

      std::vector<std::vector<std::size_t>> splitSeedSlIdxs;
      IdxTree idxTree(idxContainer, it, rootEndIt);
      constructTracks(idxTree.m_root, trackContainer, splitSeedSlIdxs);
      for (const auto& seedIdxs : splitSeedSlIdxs) {
        if (seedIdxs.size() < m_cfg.minLayers ||
            seedIdxs.size() > m_cfg.maxLayers) {
          continue;
        }

        Acts::Vector3 newDir = dir;
        Acts::Vector3 newPoint = point;
        Acts::ActsSquareMatrix<6> newCov;

        double chi2 = 0;
        for (std::size_t l = 0; l < m_cfg.nGX2Iterations; l++) {
          chi2 = m_cfg.gx2Fitter->gx2Fit(ctx.geoContext, inputSourceLinks,
                                         seedIdxs, newPoint, newDir, newCov);
        }
        if (chi2 > m_cfg.maxChi2) {
          continue;
        }

        // Transport parameteres to the reference surface
        double dVertex = (refSurfCenter - newPoint).dot(refSurfNormal) /
                         newDir.dot(refSurfNormal);
        Acts::Vector3 vertex3 = newPoint + newDir * dVertex;
        Acts::Vector4 vertex(vertex3.x(), vertex3.y(), vertex3.z(), 0);

        double thetaY = std::atan(newDir.y() / newDir.x());
        double absMom =
            std::abs(m_dipoleFieldStrength * m_dipoleLength / thetaY);

        // Estimate the abs momentum variance
        Acts::BoundMatrix originCov = m_cfg.originCov;
        originCov(Acts::eBoundQOverP, Acts::eBoundQOverP) =
            transportCovToReference(ctx.geoContext, vertex3, newPoint, newDir,
                                    newCov)(Acts::eBoundQOverP,
                                            Acts::eBoundQOverP);

        std::vector<std::size_t> sourceLinksIdxs = seedIdxs;
        sourceLinksIdxs.insert(sourceLinksIdxs.end(),
                               inputBpmSourceLinksIndices.begin(),
                               inputBpmSourceLinksIndices.end());

        if (m_cfg.propDirection == PropagationDirection::forward) {
          outSeeds.emplace_back(std::move(sourceLinksIdxs),
                                outTrackParameters.size(), outSeeds.size());
          outTrackParameters.emplace_back(vertex, newDir, 1_e / absMom,
                                          originCov,
                                          Acts::ParticleHypothesis::electron());
        } else if (m_cfg.propDirection == PropagationDirection::backward) {
          outSeeds.emplace_back(std::move(sourceLinksIdxs),
                                outTrackParameters.size(), outSeeds.size());
          outTrackParameters.emplace_back(vertex, -newDir, -1_e / absMom,
                                          originCov,
                                          Acts::ParticleHypothesis::electron());
        }
      }
    }
  }
  outSeeds.shrink_to_fit();
  outTrackParameters.shrink_to_fit();

  ACTS_DEBUG("Sending " << outSeeds.size() << " seeds");
  ACTS_DEBUG("Sending " << outTrackParameters.size() << " track parameteres");
  m_outputSeeds(ctx, std::move(outSeeds));
  m_outputTrackParameters(ctx, std::move(outTrackParameters));

  return ProcessCode::SUCCESS;
}
