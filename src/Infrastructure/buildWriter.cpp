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

  for (const auto& type : types) {
    if (!root.contains(type)) {
      throw std::runtime_error(
          "buildWriters: section ['" + type + "'] not found in config");
    }

    const auto& section = toml::find(root, type);
    writers.push_back(
        WriterRegistry::instance().build(type, section, logLevel, runRoot));
  }

  return writers;
}
}  // namespace TrackingPipeline
