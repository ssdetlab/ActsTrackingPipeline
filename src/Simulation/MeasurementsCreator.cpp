#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"
#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreatorAction.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <Acts/Utilities/VectorHelpers.hpp>
#include <cstddef>
#include <cstdlib>
#include <utility>
#include <vector>

#include <chrono>
#include <stdexcept>

MeasurementsCreator::MeasurementsCreator(
    const Propagator propagator,
    const Config& config) 
        : m_cfg(config),
        m_propagator(propagator) {};

std::tuple<std::vector<Acts::SourceLink>, SimClusters>
MeasurementsCreator::gen(
    const AlgorithmContext& ctx,
    RandomEngine& rng,
    std::size_t id) const {
        using Actions = Acts::ActionList<MeasurementsCreatorAction>;
        using Aborters = Acts::AbortList<Acts::EndOfWorldReached>;
        using PropagatorOptions =
            typename Propagator::template Options<Actions, Aborters>;

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

        // Set options for propagator
        PropagatorOptions options(
            ctx.geoContext, ctx.magFieldContext);

        auto& creator = 
            options.actionList.template get<MeasurementsCreatorAction>();
        options.maxSteps = m_cfg.maxSteps;
        creator.sourceId = id;

        // Generate initial track parameters
        rng.seed(std::chrono::system_clock::now().time_since_epoch().count());

        Acts::Vector3 spatial = m_cfg.vertexGenerator->genVertex(rng);
        Acts::Vector4 mPos4 = {spatial.x(), spatial.y(), spatial.z(), 0};

        Acts::Vector3 mom = m_cfg.momentumGenerator->genMomentum(rng);
        Acts::ActsScalar p = mom.norm(); 
        Acts::ActsScalar phi = Acts::VectorHelpers::phi(mom);
        Acts::ActsScalar theta = Acts::VectorHelpers::theta(mom);

        TrackParameters trackParameters(
            mPos4, phi, theta,
            1_e / p, ipCov, 
            Acts::ParticleHypothesis::electron());

        // Launch propagation and collect the measurements
        SimClusters simClusters;
        std::vector<Acts::SourceLink> sourceLinks;

        MeasurementsCreatorAction::result_type resultParameters;
        try {
            auto result = m_propagator.propagate(trackParameters, options).value();

            resultParameters = 
                result.template get<MeasurementsCreatorAction::result_type>();
        }
        catch (const std::runtime_error& err) {
            // std::cout << err.what() << "\n";
        }

        int trackId = (m_cfg.isSignal) ? 1 : -1;
        for (const auto& boundPars : resultParameters) {
            Acts::BoundVector boundVec = boundPars.parameters();

            Acts::GeometryIdentifier geoId = 
                boundPars.referenceSurface().geometryId();

            // Sometimes scattering makes particles 
            // to be stuck in the same surface
            if (std::ranges::find_if(simClusters,
                [&](const auto& cl) {
                    return (cl.sourceLink.geometryId() == geoId);
                }) != simClusters.end()) {
                    continue;
            }

            // Digitize hits
            Acts::Vector2 trueLocalPos{
                boundVec[Acts::eBoundLoc0], boundVec[Acts::eBoundLoc1]};

            auto [digCov, digLocalPos] = 
                m_cfg.hitDigitizer->genCluster(rng, geoId, trueLocalPos);

            // Truth information
            SimpleSourceLink hitSimpleSl(
                trueLocalPos, 
                digCov, 
                geoId, 
                ctx.eventNumber,
                -1);
            Acts::SourceLink hitSl(hitSimpleSl);

            SimHit sm{
                hitSl,
                boundVec,
                trackParameters,
                trackId,
                static_cast<int32_t>(id),
                static_cast<int32_t>(ctx.eventNumber)};

            // Observable information
            SimpleSourceLink simpleSl(
                digLocalPos, 
                digCov, 
                geoId, 
                ctx.eventNumber,
                -1);
            Acts::SourceLink sl(simpleSl);
            sourceLinks.push_back(sl);

            SimCluster cl{
                simpleSl,
                {sm},
                m_cfg.isSignal,
            };
            simClusters.push_back(cl);
        }
        
        return {std::move(sourceLinks), std::move(simClusters)};
};
