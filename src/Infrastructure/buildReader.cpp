#include "TrackingPipeline/Infrastructure/ReaderRegistry.hpp"

#include <stdexcept>
#include <string>
#include <toml.hpp>

namespace TrackingPipeline {

ReaderPtr buildReader(const std::string& type,
                      const toml::value& root,
                      const SurfaceMap& surfaceMap,
                      Acts::Logging::Level logLevel) {
  if (!root.contains(type)) {
    throw std::runtime_error(
        "buildReader: section ['" + type + "'] not found in config");
  }

  const auto& section = toml::find(root, type);

  return ReaderRegistry::instance().build(type, section, surfaceMap, logLevel);
}

}  // namespace TrackingPipeline
