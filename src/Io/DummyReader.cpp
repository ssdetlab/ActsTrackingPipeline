#include "TrackingPipeline/Io/DummyReader.hpp"

DummyReader::DummyReader(const Config& config) : IReader(), m_cfg(config) {
  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
  m_outputSimClusters.initialize(m_cfg.outputSimClusters);
  m_outputSourceLinkIndices.initialize(m_cfg.outputSourceLinkIndices);
}

ProcessCode DummyReader::read(const AlgorithmContext& ctx) {
  std::vector<Acts::SourceLink> sourceLinks{};
  std::vector<std::size_t> sourceLinksIndices{};
  SimClusters clusters{};

  m_outputSourceLinks(ctx, std::move(sourceLinks));
  m_outputSimClusters(ctx, std::move(clusters));
  m_outputSourceLinkIndices(ctx, std::move(sourceLinksIndices));

  return ProcessCode::SUCCESS;
}

std::pair<std::size_t, std::size_t> DummyReader::availableEvents() const {
  return {0, m_cfg.nEvents};
}
