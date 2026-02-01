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

#include <omp.h>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"

using go = E320Geometry::GeometryOptions;

using namespace Acts::UnitLiterals;

E320SeedingAlgorithm::E320SeedingAlgorithm(const Config& config,
                                           Acts::Logging::Level level)
    : IAlgorithm("E320SeedingAlgorithm", level), m_cfg(config) {
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_outputSeeds.initialize(m_cfg.outputSeeds);

  const auto& goInst = *go::instance();

  m_dipoleLength = goInst.dipoleHalfPrimary * 2;
  m_dipoleFieldStrength = goInst.dipoleFieldStrength;
}

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
  std::vector<std::reference_wrapper<const Acts::SourceLink>> sourceLinkRefs;
  sourceLinkRefs.reserve(inputSourceLinks.size());
  for (const auto& sl : inputSourceLinks) {
    sourceLinkRefs.push_back(std::cref(sl));
  }
  sourceLinkRefs.shrink_to_fit();

  std::vector<HoughTransformSeeder::HTSeed> htSeeds =
      m_cfg.htSeeder->findSeeds(sourceLinkRefs, m_cfg.htOptions);

  const Acts::Vector3& refSurfCenter =
      m_cfg.referenceSurface->center(ctx.geoContext);
  const Acts::Vector3& refSurfNormal = m_cfg.referenceSurface->normal(
      ctx.geoContext, refSurfCenter, Acts::Vector3::UnitX());

  ACTS_DEBUG("Found " << htSeeds.size() << " HT seeds");
  outSeeds.reserve(htSeeds.size());
  for (std::size_t i = 0; i < htSeeds.size(); i++) {
    const auto& [point, dir, cov, sl] = htSeeds.at(i);

    if (std::atan(dir.y() / dir.x() < 0.0)) {
      continue;
    }

    // Transport parameteres to the reference surface
    double dVertex =
        (refSurfCenter - point).dot(refSurfNormal) / dir.dot(refSurfNormal);
    Acts::Vector3 vertex3 = point + dir * dVertex;
    Acts::Vector4 vertex(vertex3.x(), vertex3.y(), vertex3.z(), 0);

    double thetaY = std::atan(dir.y() / dir.x());
    double absMom = std::abs(m_dipoleFieldStrength * m_dipoleLength /
                             std::sin(thetaY - m_cfg.beamlineTilt));

    outSeeds.emplace_back(sl,
                          Acts::CurvilinearTrackParameters(
                              vertex, dir, -1_e / absMom, m_cfg.originCov,
                              Acts::ParticleHypothesis::electron()),
                          i);
  }
  ACTS_DEBUG("Found " << outSeeds.size() << " seeds");
  ACTS_DEBUG("Sending " << outSeeds.size() << " seeds");
  m_outputSeeds(ctx, std::move(outSeeds));
  return ProcessCode::SUCCESS;
}
