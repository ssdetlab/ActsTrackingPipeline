#include "TrackingPipeline/Io/DummyReader.hpp"
#include "TrackingPipeline/Infrastructure/ReaderRegistry.hpp"

#include <toml.hpp>

namespace {

struct DummyReaderRegistrar {
  DummyReaderRegistrar() {
    using namespace TrackingPipeline;

    ReaderRegistry::instance().registerBuilder(
      "DummyReader",
      [](const toml::value& section,
         const SurfaceMap& /*surfaceMap*/,
         Acts::Logging::Level /*logLevel*/) -> ReaderPtr {

        DummyReader::Config cfg;

        cfg.nEvents = static_cast<std::size_t>(
          toml::find_or<int>(section, "nEvents", 1000));
        cfg.outputSourceLinks = 
            toml::find_or<std::string>(section, "outputSourceLinks", "DummySourceLinks");
        cfg.outputSimClusters = 
            toml::find_or<std::string>(section, "outputSimClusters", "DummySimClusters");

        return std::make_shared<DummyReader>(cfg);
      });
  }
} _DummyReaderRegistrar;

}  // namespace
