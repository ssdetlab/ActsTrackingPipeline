#include "TrackingPipeline/TrackFinding/TrackLookupValidationAlgorithm.hpp"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/VectorHelpers.hpp>

#include <optional>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

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

  ACTS_DEBUG("Received " << clusters.size() << " clusters");
  std::vector<Acts::CurvilinearTrackParameters> ipPars;
  std::vector<Acts::CurvilinearTrackParameters> refLayerPars;
  std::vector<Acts::CurvilinearTrackParameters> ipParsEst;
  std::vector<Acts::CurvilinearTrackParameters> refLayerParsEst;
  for (const auto& cluster : clusters) {
    ACTS_VERBOSE("Analysing cluster at GeoID "
                 << cluster.sourceLink.geometryId());
    if (!m_cfg.refLayers.contains(cluster.sourceLink.geometryId())) {
      ACTS_VERBOSE("Cluster is not in the reference layer. Continue");
      continue;
    }
    if (!cluster.isSignal) {
      ACTS_VERBOSE("Cluster is not signal. Continue");
      continue;
    }

    const Acts::Surface* layer =
        m_cfg.refLayers.at(cluster.sourceLink.geometryId());
    ACTS_VERBOSE("Cluster layer GeoID "
                 << layer->geometryId() << " at center ["
                 << layer->center(ctx.geoContext).transpose() << "]");
    for (const auto& hit : cluster.truthHits) {
      if (hit.trackId != 1) {
        ACTS_VERBOSE("Cluster hit does not belong to a signal track. Continue");
        continue;
      }

      auto ip = hit.ipParameters;

      Acts::Vector2 refLocalPos{hit.truthParameters[Acts::eBoundLoc0],
                                hit.truthParameters[Acts::eBoundLoc1]};
      Acts::Vector3 refGlobalPos = layer->localToGlobal(
          ctx.geoContext, refLocalPos, Acts::Vector3(0, 1, 0));

      Acts::Vector4 refFourPosition{refGlobalPos.x(), refGlobalPos.y(),
                                    refGlobalPos.z(), 0};

      double refPhi = hit.truthParameters[Acts::eBoundPhi];
      double refTheta = hit.truthParameters[Acts::eBoundTheta];

      Acts::Vector3 refDirection{std::sin(refTheta) * std::cos(refPhi),
                                 std::sin(refTheta) * std::sin(refPhi),
                                 std::cos(refTheta)};

      Acts::CurvilinearTrackParameters ref(
          refFourPosition, refDirection,
          hit.truthParameters[Acts::eBoundQOverP], std::nullopt,
          ip.particleHypothesis());

      SimpleSourceLink hitSsl(
          refLocalPos, refGlobalPos, cluster.sourceLink.covariance(),
          cluster.sourceLink.geometryId(), cluster.sourceLink.eventId(),
          cluster.sourceLink.index());
      auto [ipEst, refEst] =
          m_cfg.estimator(ctx.geoContext, Acts::SourceLink{hitSsl});

      ipPars.push_back(ip);
      refLayerPars.push_back(ref);
      ipParsEst.push_back(ipEst);
      refLayerParsEst.push_back(refEst);
    }
  }

  ACTS_DEBUG("Collected true " << ipPars.size() << " IP parameters");
  ACTS_DEBUG("Collected true " << refLayerPars.size()
                               << " reference layer parameters");
  ACTS_DEBUG("Estimated " << ipParsEst.size() << " IP parameters");
  ACTS_DEBUG("Estimated " << refLayerParsEst.size()
                          << " reference layer parameters");

  m_outputIpPars(ctx, std::move(ipPars));
  m_outputRefLayerPars(ctx, std::move(refLayerPars));
  m_outputIpParsEst(ctx, std::move(ipParsEst));
  m_outputRefLayerParsEst(ctx, std::move(refLayerParsEst));

  return ProcessCode::SUCCESS;
}
