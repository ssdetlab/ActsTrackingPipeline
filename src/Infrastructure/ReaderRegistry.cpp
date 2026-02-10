#include "TrackingPipeline/Infrastructure/ReaderRegistry.hpp"

#include <stdexcept>
#include <utility>

namespace TrackingPipeline {

ReaderRegistry& ReaderRegistry::instance() {
  static ReaderRegistry registry;
  return registry;
}

void ReaderRegistry::registerBuilder(const std::string& type,
                                     ReaderBuilder builder) {
  m_builders[type] = std::move(builder);
}

ReaderPtr ReaderRegistry::build(const std::string& type,
                                const toml::value& section,
                                const SurfaceMap& surfaceMap,
                                Acts::Logging::Level logLevel) const {
  auto it = m_builders.find(type);
  if (it == m_builders.end()) {
    throw std::runtime_error(
        "ReaderRegistry::build: unknown reader type '" + type + "'");
  }
  return it->second(section, surfaceMap, logLevel);
}

}  // namespace TrackingPipeline
