#include "ActsLUXEPipeline/MeasurementsCreator.hpp"

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
    const TrackParameters& trackParameters) const {
        using Actions = Acts::ActionList<MeasurementsCreatorAction>;
        using Aborters = Acts::AbortList<Acts::EndOfWorldReached>;

        // Set options for propagator
        Acts::PropagatorOptions<Actions, Aborters> options(
            ctx.geoContext, ctx.magFieldContext);

        auto& creator = options.actionList.get<MeasurementsCreatorAction>();
        creator.sourceId = ctx.eventNumber;

        // Launch and collect the measurements
        auto result = propagator.propagate(trackParameters, options).value();
        return std::move(result.template get<SimMeasurements>());
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
        Acts::Vector3 mom = m_cfg.momentumGenerator->gen(rng);

        Acts::ActsScalar p = mom.norm(); 
        Acts::ActsScalar E = std::hypot(p, me);
        Acts::ActsScalar theta = std::acos(mom.z()/p);
        Acts::ActsScalar phi = std::atan2(mom.y(), mom.x());
        
        Acts::Vector3 spatial = m_cfg.vertexGenerator->gen(rng);
        Acts::Vector4 mPos4 = {spatial.x(), spatial.y(), spatial.z(), 0};
    
        SimMeasurements results = createMeasurements(
            m_propagator, ctx,
            TrackParameters(
                mPos4, phi*180/M_PI*1_degree, 
                theta*180/M_PI*1_degree,
                -1_e / (p*1_GeV), ipCov, 
                Acts::ParticleHypothesis::electron()));
    
        m_outputMeasurements(ctx, std::move(results));
        return ProcessCode::SUCCESS;
}
