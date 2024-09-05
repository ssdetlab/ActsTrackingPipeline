#include "ActsLUXEPipeline/MeasurementsCreator.hpp"

#include <chrono>

/// @brief Constructor
MeasurementsCreator::MeasurementsCreator(
    Propagator propagator, 
    const Config& config, 
    Acts::Logging::Level level)
        : IAlgorithm("MeasurementsCreator", level),
        m_cfg(config), 
        m_propagator(std::move(propagator)) {
            m_outputMeasurements.initialize(m_cfg.outputCollection);
};
        
SimMeasurements MeasurementsCreator::createMeasurements(
    const Propagator& propagator,
    const AlgorithmContext& ctx,
    const TrackParameters& trackParameters,
    std::size_t id) const {
        using Actions = Acts::ActionList<MeasurementsCreatorAction>;
        using Aborters = Acts::AbortList<Acts::EndOfWorldReached>;
        using PropagatorOptions =
            typename Propagator::template Options<Actions, Aborters>;

        // Set options for propagator
        PropagatorOptions options(
            ctx.geoContext, ctx.magFieldContext);

        auto& creator = options.actionList.template get<MeasurementsCreatorAction>();
        creator.sourceId = id;
        options.maxSteps = 1000;

        // Launch and collect the measurements
        try {
            auto result = propagator.propagate(trackParameters, options).value();

            auto measurements = result.template get<SimMeasurements>();
    
            // Fill the true IP parameters
            for (auto& meas : measurements) {
                meas.ipParameters = trackParameters;
            }
            return std::move(measurements);
        }
        catch (const std::runtime_error& e) {
            return {};
        }
};
    
/// @brief The execute method
ProcessCode MeasurementsCreator::execute(
    const AlgorithmContext& ctx) const {
        using namespace Acts::UnitLiterals;

        // Create a random number generator
        RandomEngine rng =
            m_cfg.randomNumberSvc->spawnGenerator(ctx);

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
    
        // Create the measurements
        auto me = 0.511 * Acts::UnitConstants::MeV;
        SimMeasurements results;
        for (int i = 0; i < m_cfg.nTracks; i++) {
            rng.seed(std::chrono::system_clock::now().time_since_epoch().count());
            Acts::Vector3 mom = m_cfg.momentumGenerator->gen(rng);
    
            Acts::ActsScalar p = mom.norm(); 
            Acts::ActsScalar theta = std::acos(mom.z()/p);
            Acts::ActsScalar phi = std::atan2(mom.y(), mom.x());
            
            Acts::Vector3 spatial = m_cfg.vertexGenerator->gen(rng);
            Acts::Vector4 mPos4 = {spatial.x(), spatial.y(), spatial.z(), 0};

            SimMeasurements res = createMeasurements(
                m_propagator, ctx,
                TrackParameters(
                    mPos4, phi, theta,
                    1_e / p, ipCov, 
                    Acts::ParticleHypothesis::electron()),
                    i);

            results.insert(results.end(), res.begin(), res.end());
        }

        m_outputMeasurements(ctx, std::move(results));
        return ProcessCode::SUCCESS;
}
