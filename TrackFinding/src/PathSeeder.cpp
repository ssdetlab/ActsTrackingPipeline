#include "ActsLUXEPipeline/PathSeeder.hpp"

#include "ActsLUXEPipeline/SimpleSourceLink.hpp"

#include "Acts/Definitions/Units.hpp"

Acts::ActsScalar PathSeeder::findClosestValue(
    std::unordered_map<
        Acts::ActsScalar,
        Acts::ActsScalar>& lookupTable, 
    Acts::ActsScalar x) const {
        Acts::ActsScalar closestKey;
        Acts::ActsScalar closestValue;
        bool first = true;
    
        for (const auto& entry : lookupTable) {
            if (first) {
                closestKey = entry.first;
                closestValue = entry.second;
                first = false;
            }
            if (std::abs(entry.first - x) < std::abs(closestKey - x)) {
                closestKey = entry.first;
                closestValue = entry.second;
            }
        }
        return closestValue;
}

ProcessCode PathSeeder::execute(const AlgorithmContext& ctx) const {
    using namespace Acts::UnitLiterals;

    // Get the input measurements
    // from the context
    auto input = m_inputMeasurements(ctx);

    if (input.empty()) {
        m_outputSeeds(ctx, Seeds());
        return ProcessCode::SUCCESS;
    }

    // Get the input source links
    std::vector<Acts::SourceLink> inpSourceLinks;
    for (const auto& meas : input) {
        inpSourceLinks.push_back(meas.sourceLink);
    }

    // Bin the source links
    m_cfg.sourceLinkBinner->initialize(ctx.geoContext, inpSourceLinks);

    // Read lookup tables
    auto EXFirstLookUp = m_cfg.EXFirstLookupProvider->getLookup();
    auto XFirstXLastLookUp = m_cfg.XFirstXLastLookupProvider->getLookup();
    auto XFisrtYLastLookUp = m_cfg.XFirstYLastLookupProvider->getLookup();
    auto ZFirstZLastLookUp = m_cfg.ZFirstZLastLookupProvider->getLookup();

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

    // Create the seeds
    Seeds seeds;
    std::vector<Acts::SourceLink> sourceLinks;
    for (int i = 0; i < input.size() - 1; i++) {
        auto ssl = input.at(i).sourceLink.get<SimpleSourceLink>();

        Acts::Vector3 globalPos = m_cfg.surfaceAccessor(
            input.at(i).sourceLink)->localToGlobal(
                ctx.geoContext, 
                ssl.parameters, 
                Acts::Vector3{0, 1, 0});

        Acts::ActsScalar x1 = globalPos[0];
        Acts::ActsScalar y1 = globalPos[1];
        Acts::ActsScalar z1 = globalPos[2];

        // Check if the hit is in the first layer
        if (!m_cfg.firstLayerExtent.contains(globalPos)) {
            continue;
        }

        // Get the energy and the last hit position
        // from the lookup tables 
        Acts::ActsScalar E = findClosestValue(EXFirstLookUp, x1);
        Acts::ActsScalar x4 = findClosestValue(XFirstXLastLookUp, x1);
        Acts::ActsScalar z4 = findClosestValue(ZFirstZLastLookUp, z1);
        Acts::ActsScalar y4 = findClosestValue(XFisrtYLastLookUp, x1);

        Acts::ActsScalar me = 0.511 * Acts::UnitConstants::MeV;
        // Momentum magnitude
        Acts::ActsScalar pMagnitude = 
            std::sqrt(std::pow(E, 2) - std::pow(me, 2));

        // Direction vector
        auto dir = Acts::Vector3(
            (x4 - x1), (y4 - y1), (z4 - z1));
        dir.normalize();

        auto dirIp = Acts::Vector3(0, sqrt(1 - dir.z()*dir.z()), dir.z());

        // Ip parameter is the same for all hits
        // with the same track id
        Acts::CurvilinearTrackParameters ipParameters(
            Acts::Vector4(0, 0, 0, 0),
            dirIp,
            -1_e/pMagnitude,
            ipCov,
            Acts::ParticleHypothesis::electron());

        // Intersect with the surfaces
        std::vector<Acts::SurfaceIntersection> intersections =
            m_cfg.intersectionFinder->findIntersections(
                ctx.geoContext, globalPos, dir);

        if (intersections.empty()) {
            continue;
        }

        // Store the pivot source link
        sourceLinks.push_back(input.at(i).sourceLink);

        // Iterate over the intersections 
        // and get the source links
        for (auto intersect : intersections) {
            // Get the intersection point
            Acts::Vector3 refPoint = intersect.position();

            // Get the path width
            auto [pathWidthX, pathWidthZ] = 
                m_cfg.pathWidthProvider->getPathWidth(
                    ctx.geoContext, 
                    intersect.object()->geometryId());

            // Get the bounds of the path
            Acts::ActsScalar topX = 
                refPoint.x() + pathWidthX;
            Acts::ActsScalar botX = 
                refPoint.x() - pathWidthX;
            Acts::ActsScalar topZ = 
                refPoint.z() + pathWidthZ;
            Acts::ActsScalar botZ = 
                refPoint.z() - pathWidthZ;

            // Get the lookup table for the source links
            auto grid = m_cfg.sourceLinkBinner->getLookupTable(
                intersect.object()->geometryId());

            // Get the range of bins to search for source links
            auto botLeftBin = grid.localBinsFromPosition(
                Acts::Vector2(botX, botZ));
            auto topRightBin = grid.localBinsFromPosition(
                Acts::Vector2(topX, topZ));

            // Get the source links from the lookup table
            // by iterating over the bin ranges
            auto currentBin = botLeftBin;
            while (currentBin.at(1) <= topRightBin.at(1)) {
                while (currentBin.at(0) <= topRightBin.at(0)) {
                    auto sourceLinksToAdd = 
                        grid.atLocalBins(currentBin);

                    sourceLinks.insert(
                        sourceLinks.end(),
                        sourceLinksToAdd.begin(),
                        sourceLinksToAdd.end());
                    currentBin.at(0)++;
                }
                currentBin.at(1)++;
                currentBin.at(0) = botLeftBin.at(0);
            }
        }

        // Add the seed to the list
        seeds.push_back(Seed
            {sourceLinks, ipParameters, i});
        sourceLinks.clear();
    }

    m_outputSeeds(ctx, std::move(seeds));

    return ProcessCode::SUCCESS;
}
