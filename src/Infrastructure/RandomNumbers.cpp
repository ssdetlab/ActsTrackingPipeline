#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"

RandomNumbers::RandomNumbers(const Config& cfg) : m_cfg(cfg) {}

RandomEngine RandomNumbers::spawnGenerator() const {
  return RandomEngine(
      std::chrono::system_clock::now().time_since_epoch().count());
}

RandomEngine RandomNumbers::spawnGenerator(
    const AlgorithmContext& context) const {
  return RandomEngine(generateSeed(context));
}

uint64_t RandomNumbers::generateSeed(const AlgorithmContext& context) const {
  return m_cfg.seed + context.eventNumber;
}
