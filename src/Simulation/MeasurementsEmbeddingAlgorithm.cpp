#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"

#include <stdexcept>

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

  ACTS_DEBUG("Received " << clusters.size() << " clusters");
  ACTS_DEBUG("Received " << sourceLinks.size() << " source links");

  RandomEngine rng = m_cfg.randomNumberSvc->spawnGenerator();

  // Create the measurements
  ACTS_DEBUG("Starting propagation of " << m_cfg.nMeasurements << " tracks");
  std::size_t inputSize = sourceLinks.size();
  std::size_t outputSize = inputSize + m_cfg.nMeasurements;
  sourceLinks.reserve(outputSize);
  int index = inputSize;
  for (std::size_t i = inputSize; i < outputSize; i++) {
    const auto& [sls, cls] = m_cfg.measurementGenerator->gen(ctx, rng, index);
    index += sls.size();

    ACTS_VERBOSE("Created " << sls.size() << " measurements");
    sourceLinks.insert(sourceLinks.end(), sls.begin(), sls.end());
    clusters.insert(clusters.end(), cls.begin(), cls.end());
  }

  ACTS_DEBUG("Created " << sourceLinks.size() << " measurements in total");
  m_outputSourceLinks(ctx, std::move(sourceLinks));
  m_outputSimClusters(ctx, std::move(clusters));

  return ProcessCode::SUCCESS;
}
