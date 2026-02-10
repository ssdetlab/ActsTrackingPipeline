#include "TrackingPipeline/TrackFinding/PathSeedingAlgorithm.hpp"

/// @brief Constructor
PathSeedingAlgorithm::PathSeedingAlgorithm(
    const Config& config, Acts::Logging::Level level)
    : IAlgorithm("PathSeedingAlgorithm", level),
    m_cfg(config) {
        m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
        m_outputSeeds.initialize(m_cfg.outputSeeds);
}

/// @brief The execute method        
ProcessCode PathSeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
        using namespace Acts::UnitLiterals;

        // Get the input measurements
        // from the context
        auto input = m_inputSourceLinks(ctx);

        if (input.empty()) {
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
        Acts::BoundSquareMatrix ipCov =
            ipStdDev.cwiseProduct(ipStdDev).asDiagonal();

        auto gridLookup = m_cfg.sourceLinkGridConstructor->constructGrid(
            ctx.geoContext, input);

        std::vector<Acts::PathSeeder::PathSeed> pathSeeds;
        m_cfg.seeder->findSeeds(
            ctx.geoContext, 
            gridLookup, 
            pathSeeds);

        Seeds outSeeds;
        auto me = 0.511 * Acts::UnitConstants::MeV;

        for (std::int32_t i = 0; i < pathSeeds.size(); i++) {
            const auto& seed = pathSeeds.at(i);

            if (seed.second.size() < m_cfg.minSeedSize || 
                seed.second.size() > m_cfg.maxSeedSize) {
                    continue;
            }

            outSeeds.push_back(Seed{
                seed.second,
                seed.first,
                i});
        }
    
        m_outputSeeds(ctx, std::move(outSeeds));
        
        return ProcessCode::SUCCESS;
}