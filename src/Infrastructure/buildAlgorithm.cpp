#include "TrackingPipeline/Infrastructure/AlgorithmRegistry.hpp"

#include <stdexcept>
#include <string>
#include <vector>
#include <toml.hpp>

namespace TrackingPipeline {

std::vector<AlgorithmPtr> buildAlgorithms(
    const std::vector<std::string>& types,
    const toml::value& root,
    Acts::Logging::Level logLevel) {

  std::vector<AlgorithmPtr> algos;
  algos.reserve(types.size());

  for (const auto& type : types) {
    if (!root.contains(type)) {
      throw std::runtime_error(
          "buildAlgorithms: section ['" + type + "'] not found in config");
    }
    const auto& section = toml::find(root, type);
    algos.push_back(
        AlgorithmRegistry::instance().build(type, section, logLevel));
  }

  return algos;
}
}  // namespace TrackingPipeline
