#include "TrackingPipeline/Infrastructure/WriterRegistry.hpp"

#include <stdexcept>
#include <utility>

namespace TrackingPipeline {

WriterRegistry& WriterRegistry::instance() {
  static WriterRegistry registry;
  return registry;
}

void WriterRegistry::registerBuilder(const std::string& type,
                                     WriterBuilder builder) {
  m_builders[type] = std::move(builder);
}

WriterPtr WriterRegistry::build(const std::string& type,
                                const toml::value& section,
                                Acts::Logging::Level logLevel,
                                const std::string& runRoot) const {
  auto it = m_builders.find(type);
  if (it == m_builders.end()) {
    throw std::runtime_error(
        "WriterRegistry::build: unknown writer type '" + type + "'");
  }
  return it->second(section, logLevel, runRoot);
}

}  // namespace TrackingPipeline
