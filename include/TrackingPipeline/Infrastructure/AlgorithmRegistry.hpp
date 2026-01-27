#pragma once

#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <toml.hpp>

namespace TrackingPipeline {

using AlgorithmPtr = std::shared_ptr<IAlgorithm>;

using AlgorithmBuilder =
    std::function<AlgorithmPtr(const toml::value& section,
                               Acts::Logging::Level logLevel)>;
class AlgorithmRegistry {
 public:
  static AlgorithmRegistry& instance();

  void registerBuilder(const std::string& type, AlgorithmBuilder builder);

  AlgorithmPtr build(const std::string& type,
                     const toml::value& section,
                     Acts::Logging::Level logLevel) const;

 private:
  AlgorithmRegistry() = default;
  AlgorithmRegistry(const AlgorithmRegistry&) = delete;
  AlgorithmRegistry& operator=(const AlgorithmRegistry&) = delete;
  std::unordered_map<std::string, AlgorithmBuilder> m_builders;
};

std::vector<AlgorithmPtr> buildAlgorithms(
    const std::vector<std::string>& types,
    const toml::value& root,
    Acts::Logging::Level logLevel);

}  // namespace TrackingPipeline
