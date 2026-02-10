#include "TrackingPipeline/Infrastructure/WriterRegistry.hpp"

#include <stdexcept>
#include <string>
#include <vector>
#include <toml.hpp>

namespace TrackingPipeline {

std::vector<WriterPtr> buildWriters(const std::vector<std::string>& types,
                                    const toml::value& root,
                                    Acts::Logging::Level logLevel,
                                    const std::string& runRoot) {

  std::vector<WriterPtr> writers;
  writers.reserve(types.size());

  for (const auto& name : types) {
    if (!root.contains(name)) {
      throw std::runtime_error(
          "buildWriters: section ['" + name + "'] not found in config");
    }

    const auto& section = toml::find(root, name);
    if (!section.contains("type")) {
      throw std::runtime_error(
          "buildWriters: section ['" + name +
          "'] has no 'type' field (expected writer type string)");
    }
    const auto writerType = toml::find<std::string>(section, "type");
    writers.push_back(
        WriterRegistry::instance().build(writerType, section, logLevel, runRoot));
  }

  return writers;
}
}  // namespace TrackingPipeline
