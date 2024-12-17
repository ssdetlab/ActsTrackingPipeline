#include "TrackingPipeline/TrackFinding/TrackLookupValidationAlgorithm.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/VectorHelpers.hpp>
#include <optional>

TrackLookupValidationAlgorithm::TrackLookupValidationAlgorithm(
    const Config& config, Acts::Logging::Level level)
    : IAlgorithm("TrackParamsLookupValidation", level),
    m_cfg(std::move(config)) {
        m_inputClusters.initialize(m_cfg.inputClusters);

        m_outputIpPars.initialize(m_cfg.outputIpPars);
        m_outputRefLayerPars.initialize(m_cfg.outputRefLayerPars);
        m_outputIpParsEst.initialize(m_cfg.outputIpParsEst);
        m_outputRefLayerParsEst.initialize(m_cfg.outputRefLayerParsEst);
}

ProcessCode TrackLookupValidationAlgorithm::execute(
    const AlgorithmContext& ctx) const {
        using namespace Acts::UnitLiterals;

        auto clusters = m_inputClusters(ctx);

        std::vector<Acts::CurvilinearTrackParameters> ipPars;
        std::vector<Acts::CurvilinearTrackParameters> refLayerPars;
        std::vector<Acts::CurvilinearTrackParameters> ipParsEst;
        std::vector<Acts::CurvilinearTrackParameters> refLayerParsEst;
        for (const auto& cluster : clusters) {
            if (!m_cfg.refLayers.contains(cluster.sourceLink.geometryId())) {
                continue;
            }
            if (!cluster.isSignal) {
                continue;
            }

            const Acts::Surface* layer = m_cfg.refLayers.at(cluster.sourceLink.geometryId());
            for (const auto& hit : cluster.truthHits) {
                if (hit.trackId != 1) {
                    continue;
                }
                
                auto ip = hit.ipParameters;

                Acts::Vector2 refLocalPos{
                    hit.truthParameters[Acts::eBoundLoc0],
                    hit.truthParameters[Acts::eBoundLoc1]};
                Acts::Vector3 refGlobalPos = layer->localToGlobal(
                    ctx.geoContext, 
                    refLocalPos, 
                    Acts::Vector3(0, 1, 0));

                Acts::Vector4 refFourPosition{
                    refGlobalPos.x(),
                    refGlobalPos.y(),
                    refGlobalPos.z(),
                    0};

                Acts::ActsScalar refPhi = hit.truthParameters[Acts::eBoundPhi];
                Acts::ActsScalar refTheta = hit.truthParameters[Acts::eBoundTheta];

                Acts::Vector3 refDirection{
                    std::sin(refTheta) * std::cos(refPhi),
                    std::sin(refTheta) * std::sin(refPhi),                    
                    std::cos(refTheta)};
                
                Acts::CurvilinearTrackParameters ref(
                    refFourPosition,
                    refDirection,
                    hit.truthParameters[Acts::eBoundQOverP],
                    std::nullopt,
                    ip.particleHypothesis());
                
                auto [ipEst, refEst] =
                    m_cfg.estimator(ctx.geoContext, hit.sourceLink);

                ipPars.push_back(ip);
                refLayerPars.push_back(ref);
                ipParsEst.push_back(ipEst);
                refLayerParsEst.push_back(refEst);
            }
        }
        
        m_outputIpPars(ctx, std::move(ipPars));
        m_outputRefLayerPars(ctx, std::move(refLayerPars));
        m_outputIpParsEst(ctx, std::move(ipParsEst));
        m_outputRefLayerParsEst(ctx, std::move(refLayerParsEst));
        
        return ProcessCode::SUCCESS;
}
