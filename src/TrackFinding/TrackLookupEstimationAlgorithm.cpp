#include "TrackingPipeline/TrackFinding/TrackLookupEstimationAlgorithm.hpp"

#include "Acts/Surfaces/RectangleBounds.hpp"

#include "TrackingPipeline/Infrastructure/AlgorithmRegistry.hpp"
#include "TrackingPipeline/Io/JsonTrackLookupWriter.hpp"
#include "TrackingPipeline/TrackFitting/FittingServices.hpp"
#include <toml.hpp>

TrackLookupEstimationAlgorithm::TrackLookupEstimationAlgorithm(
    const Config& config, Acts::Logging::Level level)
    : IAlgorithm("TrackLookupEstimationAlgorithm", level), m_cfg(config) {
  // Iterate over the reference layers and create
  // track parameter accumulators
  for (const auto& [geoId, refSurface] : m_cfg.refLayers) {
    // Get bounds to construct the accumulator grid
    auto bounds =
        dynamic_cast<const Acts::RectangleBounds*>(&refSurface->bounds());

    if (bounds == nullptr) {
      throw std::invalid_argument("Only rectangle bounds supported");
    }
    if (refSurface->type() != Acts::Surface::SurfaceType::Plane) {
      throw std::invalid_argument("Only plane surfaces supported");
    }

    // Initialize the accumulator grid
    auto halfX = bounds->halfLengthX();
    auto halfY = bounds->halfLengthY();

    TrackLookupAxisGen axisGen{
        {-halfX, halfX}, m_cfg.bins.first, {-halfY, halfY}, m_cfg.bins.second};

    // Each reference layer has its own accumulator
    m_accumulators[geoId] =
        std::make_unique<TrackLookupAccumulator>(TrackLookupGrid(axisGen()));
  }

  m_inputClusters.initialize(m_cfg.inputClusters);
}

ProcessCode TrackLookupEstimationAlgorithm::finalize() {
  // Finiliaze the lookup tables and write them
  TrackLookup lookup;
  for (auto& [id, acc] : m_accumulators) {
    lookup.insert({id, acc->finalizeLookup()});
  }
  for (const auto& writer : m_cfg.trackLookupGridWriters) {
    writer->writeLookup(lookup);
  }

  return ProcessCode::SUCCESS;
};

ProcessCode TrackLookupEstimationAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Get the particles and hits
  const auto& clusters = m_inputClusters(ctx);

  // Iterate over the reference layer hits and
  // accumulate the track parameters
  for (const auto& cluster : clusters) {
    const auto& geoId = cluster.sourceLink.geometryId();
    if (!m_cfg.refLayers.contains(geoId)) {
      continue;
    }
    const auto& layer = m_cfg.refLayers.at(geoId);
    for (const auto& hit : cluster.truthHits) {
      const auto& ip = hit.ipParameters;

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

      // Add the track parameters to the accumulator grid
      m_accumulators.at(geoId)->addTrack(ip, ref, refLocalPos);
    }
  }

  return ProcessCode::SUCCESS;
}


namespace {

struct TrackLookupEstimationAlgorithmRegistrar {
  TrackLookupEstimationAlgorithmRegistrar() {
    using namespace TrackingPipeline;

    AlgorithmRegistry::instance().registerBuilder(
        "TrackLookupEstimationAlgorithm",
        [](const toml::value& section,
           Acts::Logging::Level logLevel) -> AlgorithmPtr {
          TrackLookupEstimationAlgorithm::Config cfg;

          // 1) Ref layers from FittingServices
          auto& fs = FittingServices::instance();
          cfg.refLayers = fs.lookupRefLayers;

          // 2) Binning from TOML
          const auto binsX = toml::find<std::size_t>(section, "binsX");
          const auto binsY = toml::find<std::size_t>(section, "binsY");
          cfg.bins = {binsX, binsY};

          // 3) Input clusters name 
          cfg.inputClusters =
              toml::find<std::string>(section, "inputClusters");

          // 4) Build the JsonTrackLookupWriter helper
          JsonTrackLookupWriter::Config writerCfg;
          writerCfg.path =
              toml::find<std::string>(section, "lookupJsonPath");

          auto jsonWriter =
              std::make_shared<JsonTrackLookupWriter>(writerCfg);

          cfg.trackLookupGridWriters = {jsonWriter};

          return std::make_shared<TrackLookupEstimationAlgorithm>(cfg, logLevel);
        });
  }
} _TrackLookupEstimationAlgorithmRegistrar;

}  // namespace