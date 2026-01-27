#include "TrackingPipeline/Infrastructure/AlgorithmRegistry.hpp"

#include <stdexcept>
#include <utility>

namespace TrackingPipeline {

AlgorithmRegistry& AlgorithmRegistry::instance() {
  static AlgorithmRegistry reg;
  return reg;
}

void AlgorithmRegistry::registerBuilder(const std::string& type,
                                      AlgorithmBuilder builder) {
  m_builders[type] = std::move(builder);
}

AlgorithmPtr AlgorithmRegistry::build(const std::string& type,
                                      const toml::value& section,
                                      Acts::Logging::Level logLevel) const {
  auto it = m_builders.find(type);
  if (it == m_builders.end()) {
    throw std::runtime_error(
        "AlgorithmRegistry::build: unknown algorithm type '" + type + "'");
  }
  return it->second(section, logLevel);
}

}  // namespace TrackingPipeline
