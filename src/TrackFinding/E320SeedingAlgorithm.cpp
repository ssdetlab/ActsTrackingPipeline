#include "TrackingPipeline/TrackFinding/E320SeedingAlgorithm.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>

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

  Acts::BoundVector ipStdDev;
  ipStdDev[Acts::eBoundLoc0] = 100_um;
  ipStdDev[Acts::eBoundLoc1] = 100_um;
  ipStdDev[Acts::eBoundTime] = 25_ns;
  ipStdDev[Acts::eBoundPhi] = 2_degree;
  ipStdDev[Acts::eBoundTheta] = 2_degree;
  ipStdDev[Acts::eBoundQOverP] = 1 / 100_GeV;
  m_ipCov = ipStdDev.cwiseProduct(ipStdDev).asDiagonal();

  const auto& goInst = *go::instance();

  m_dipoleLength = goInst.dipoleHalfPrimary * 2;
  m_dipoleFieldStrength = goInst.dipoleFieldStrength;

  m_xCorrectorLength = goInst.xCorrectorHalfPrimary * 2;
  m_xCorrectorFieldStrength = goInst.xCorrectorFieldStrength;

  m_detFirstLayerPoint[goInst.primaryIdx] = goInst.ipTcDistance;
  m_detFirstLayerPoint[goInst.longIdx] = goInst.tcCenterLong;
  m_detFirstLayerPoint[goInst.shortIdx] = goInst.tcCenterShort;
  m_detFirstLayerNormal = goInst.primaryDir;

  m_backShift = Acts::Vector3::Zero();
  m_backShift[goInst.primaryIdx] = 0.3;
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

  ACTS_DEBUG("Found " << htSeeds.size() << " HT seeds");
  outSeeds.reserve(htSeeds.size());
#pragma omp parallel for num_threads(32)
  for (std::size_t i = 0; i < htSeeds.size(); i++) {
    const auto& [point, dir, sl] = htSeeds.at(i);

    double dVertex = (m_detFirstLayerPoint - m_backShift - point)
                         .dot(m_detFirstLayerNormal) /
                     dir.dot(m_detFirstLayerNormal);

    Acts::Vector3 vertex3 = point + dir * dVertex;
    Acts::Vector4 vertex(vertex3.x(), vertex3.y(), vertex3.z(), 0);

    double thetaY = std::atan(dir.y() / dir.x());
    double thetaZ = std::atan(dir.z() / dir.x());

    double momFromDipole =
        std::abs(m_dipoleFieldStrength * m_dipoleLength / std::sin(thetaY));

    double absMom = momFromDipole;

#pragma omp critical
    {
      outSeeds.emplace_back(sl,
                            Acts::CurvilinearTrackParameters(
                                vertex, dir, -1_e / absMom, m_ipCov,
                                Acts::ParticleHypothesis::electron()),
                            i);
    }
  }
  ACTS_DEBUG("Found " << outSeeds.size() << " seeds");
  ACTS_DEBUG("Sending " << outSeeds.size() << " seeds");
  m_outputSeeds(ctx, std::move(outSeeds));
  return ProcessCode::SUCCESS;
}
