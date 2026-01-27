#pragma once

#include "TrackingPipeline/Infrastructure/IWriter.hpp"

#include "Acts/Utilities/Logger.hpp"

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

// TOML
#include <toml.hpp>

namespace TrackingPipeline {

using WriterPtr = std::shared_ptr<IWriter>;

// Builder takes a TOML section and a log level and returns a constructed writer.
using WriterBuilder =
    std::function<WriterPtr(const toml::value& section,
                            Acts::Logging::Level logLevel,
                            const std::string& runRoot)>;

class WriterRegistry {
 public:
  static WriterRegistry& instance();

  void registerBuilder(const std::string& type, WriterBuilder builder);

  WriterPtr build(const std::string& type,
                  const toml::value& section,
                  Acts::Logging::Level logLevel,
                  const std::string& runRoot) const;

 private:
  WriterRegistry() = default;
  WriterRegistry(const WriterRegistry&) = delete;
  WriterRegistry& operator=(const WriterRegistry&) = delete;

  std::unordered_map<std::string, WriterBuilder> m_builders;
};

std::vector<WriterPtr> buildWriters(const std::vector<std::string>& types,
                                    const toml::value& root,
                                    Acts::Logging::Level logLevel,
                                    const std::string& runRoot);

}  // namespace TrackingPipeline
