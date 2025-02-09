#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"

#include "Acts/EventData/SourceLink.hpp"

#include <stdexcept>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"

MeasurementsEmbeddingAlgorithm::MeasurementsEmbeddingAlgorithm(
    const Config& config, Acts::Logging::Level level)
    : IAlgorithm("MeasurementsEmbeddingAlgorithm", level), m_cfg(config) {
  if (m_cfg.measurementGenerator == nullptr) {
    throw std::runtime_error("Generator is not initialized");
  }

  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_inputSimClusters.initialize(m_cfg.inputSimClusters);

  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
  m_outputSimClusters.initialize(m_cfg.outputSimClusters);
}

ProcessCode MeasurementsEmbeddingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Get the inputs from the context
  auto sourceLinks = m_inputSourceLinks(ctx);
  auto clusters = m_inputSimClusters(ctx);

  // Create a random number generator
  RandomEngine rng = m_cfg.randomNumberSvc->spawnGenerator(ctx);

  // Create the measurements
  for (int i = 0; i < m_cfg.nMeasurements; i++) {
    auto [sls, cls] = m_cfg.measurementGenerator->gen(ctx, rng, i);

    for (std::size_t i = 0; i < sls.size(); i++) {
      auto ssl = sls.at(i).get<SimpleSourceLink>();
      ssl.setIndex(sourceLinks.size());

      auto cl = cls.at(i);
      cl.sourceLink.setIndex(ssl.index());

      if (m_cfg.clusterFilter == nullptr ||
          m_cfg.clusterFilter->operator()(ctx.geoContext, cls.at(i))) {
        sourceLinks.push_back(Acts::SourceLink(ssl));
        clusters.push_back(cl);
      }
    }
  }

  m_outputSourceLinks(ctx, std::move(sourceLinks));
  m_outputSimClusters(ctx, std::move(clusters));

  return ProcessCode::SUCCESS;
}
