#include "TrackingPipeline/TrackFinding/PathSeedingAlgorithm.hpp"

#include <Acts/Utilities/Logger.hpp>

#include <cstddef>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

/// @brief Constructor
PathSeedingAlgorithm::PathSeedingAlgorithm(const Config& config,
                                           Acts::Logging::Level level)
    : IAlgorithm("PathSeedingAlgorithm", level), m_cfg(config) {
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_outputSeeds.initialize(m_cfg.outputSeeds);
}

/// @brief The execute method
ProcessCode PathSeedingAlgorithm::execute(const AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

  // Get the input measurements
  // from the context
  auto input = m_inputSourceLinks(ctx);

  ACTS_DEBUG("Received " << input.size() << " source links");

  if (input.empty()) {
    ACTS_DEBUG("Input is empty. Skipping");
    m_outputSeeds(ctx, Seeds());
    return ProcessCode::SUCCESS;
  }

  // Create IP covariance matrix from
  // reasonable standard deviations
  Acts::BoundVector ipStdDev;
  ipStdDev[Acts::eBoundLoc0] = 100_um;
  ipStdDev[Acts::eBoundLoc1] = 100_um;
  ipStdDev[Acts::eBoundTime] = 25_ns;
  ipStdDev[Acts::eBoundPhi] = 2_degree;
  ipStdDev[Acts::eBoundTheta] = 2_degree;
  ipStdDev[Acts::eBoundQOverP] = 1 / 100_GeV;
  Acts::BoundSquareMatrix ipCov = ipStdDev.cwiseProduct(ipStdDev).asDiagonal();

  auto gridLookup =
      m_cfg.sourceLinkGridConstructor->constructGrid(ctx.geoContext, input);

  std::vector<Acts::PathSeeder::PathSeed> pathSeeds;
  m_cfg.seeder->findSeeds(ctx.geoContext, gridLookup, pathSeeds);
  ACTS_DEBUG("Found " << pathSeeds.size() << " seeds");

  Seeds outSeeds;
  auto me = 0.511 * Acts::UnitConstants::MeV;

  for (std::size_t i = 0; i < pathSeeds.size(); i++) {
    const auto& seed = pathSeeds.at(i);

    std::set<Acts::GeometryIdentifier> seedGeoIds;
    for (const auto& sl : seed.second) {
      seedGeoIds.insert(sl.get<SimpleSourceLink>().geometryId());
    }
    if (seedGeoIds.size() < m_cfg.minLayers ||
        seedGeoIds.size() > m_cfg.maxLayers) {
      continue;
    }
    if (seed.second.size() < m_cfg.minSeedSize ||
        seed.second.size() > m_cfg.maxSeedSize) {
      continue;
    }

    outSeeds.push_back(Seed{seed.second, seed.first, static_cast<int>(i)});
  }
  ACTS_DEBUG("Sending " << outSeeds.size() << " seeds");
  m_outputSeeds(ctx, std::move(outSeeds));

  return ProcessCode::SUCCESS;
}
