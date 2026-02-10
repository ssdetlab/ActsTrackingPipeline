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

  for (const auto& name : types) {
    if (!root.contains(name)) {
      throw std::runtime_error(
          "buildAlgorithms: section ['" + name + "'] not found in config");
    }
    const auto& section = toml::find(root, name);

    if (!section.contains("type")) {
      throw std::runtime_error(
          "buildAlgorithms: section ['" + name +
          "'] has no 'type' field (expected algorithm type string)");
    }
    const auto algoType = toml::find<std::string>(section, "type");
    algos.push_back(
        AlgorithmRegistry::instance().build(algoType, section, logLevel));
  }

  return algos;
}
}  // namespace TrackingPipeline
