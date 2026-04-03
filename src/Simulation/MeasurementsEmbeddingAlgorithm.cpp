#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"

#include <cstddef>
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
  m_inputSourceLinksIndices.initialize(m_cfg.inputSourceLinkIndices);

  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
  m_outputSimClusters.initialize(m_cfg.outputSimClusters);
  m_outputSourceLinkIndices.initialize(m_cfg.outputSourceLinkIndices);
}

ProcessCode MeasurementsEmbeddingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Get the inputs from the context
  auto inputSourceLinks = m_inputSourceLinks(ctx);
  auto inputSimClusters = m_inputSimClusters(ctx);
  auto inputSourceLinksIndices = m_inputSourceLinksIndices(ctx);

  ACTS_DEBUG("Received " << inputSimClusters.size() << " clusters");
  ACTS_DEBUG("Received " << inputSourceLinks.size() << " source links");
  ACTS_DEBUG("Received " << inputSourceLinksIndices.size()
                         << " source link indices");

  RandomEngine rng = m_cfg.randomNumberSvc->spawnGenerator();

  // Create the measurements
  ACTS_DEBUG("Starting propagation of " << m_cfg.nMeasurements << " tracks");
  std::size_t inputSize = inputSourceLinks.size();
  std::size_t outputSize = inputSize + m_cfg.nMeasurements;
  inputSourceLinks.reserve(outputSize);
  inputSimClusters.reserve(outputSize);
  inputSourceLinksIndices.reserve(outputSize);
  int index = inputSize;
  for (std::size_t i = inputSize; i < outputSize; i++) {
    std::size_t resSize = m_cfg.measurementGenerator->gen(
        ctx, rng, index, inputSourceLinks, inputSimClusters,
        inputSourceLinksIndices);
    index += resSize;

    ACTS_VERBOSE("Created " << resSize << " measurements");
  }

  ACTS_DEBUG("Created " << inputSourceLinks.size() << " measurements in total");
  m_outputSourceLinks(ctx, std::move(inputSourceLinks));
  m_outputSimClusters(ctx, std::move(inputSimClusters));
  m_outputSourceLinkIndices(ctx, std::move(inputSourceLinksIndices));

  return ProcessCode::SUCCESS;
}
