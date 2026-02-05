#pragma once

#include "TrackingPipeline/Infrastructure/IReader.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <toml.hpp>

using SurfaceMap = std::map<Acts::GeometryIdentifier, const Acts::Surface*>;

namespace TrackingPipeline {

using ReaderPtr = std::shared_ptr<IReader>;

using ReaderBuilder =
    std::function<ReaderPtr(const toml::value& section,
                            const SurfaceMap& surfaceMap,
                            Acts::Logging::Level logLevel)>;


class ReaderRegistry {
 public:
  static ReaderRegistry& instance();

  void registerBuilder(const std::string& type, ReaderBuilder builder);


  ReaderPtr build(const std::string& type,
                  const toml::value& section,
                  const SurfaceMap& surfaceMap,
                  Acts::Logging::Level logLevel) const;

 private:
  ReaderRegistry() = default;
  ReaderRegistry(const ReaderRegistry&) = delete;
  ReaderRegistry& operator=(const ReaderRegistry&) = delete;

  std::unordered_map<std::string, ReaderBuilder> m_builders;
};

ReaderPtr buildReader(const std::string& type,
                      const toml::value& root,
                      const SurfaceMap& surfaceMap,
                      Acts::Logging::Level logLevel);

}  // namespace TrackingPipeline
