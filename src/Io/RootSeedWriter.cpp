#include "TrackingPipeline/Io/RootSeedWriter.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <ranges>
#include <vector>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Infrastructure/WriterRegistry.hpp"

#include <toml.hpp>

RootSeedWriter::RootSeedWriter(const Config& config, Acts::Logging::Level level)
    : m_cfg(config), m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_file = new TFile(m_cfg.filePath.c_str(), "RECREATE");
  m_tree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());

  //------------------------------------------------------------------
  // Tree branches
  int buf_size = 32000;
  int split_lvl = 0;

  // Measurements
  m_tree->Branch("measurementsGlob", &m_seedMeasurementsGlob, buf_size,
                 split_lvl);
  m_tree->Branch("measurementsLoc", &m_seedMeasurementsLoc, buf_size,
                 split_lvl);
  m_tree->Branch("geoIds", &m_geoIds, buf_size, split_lvl);

  // Seed properties
  m_tree->Branch("eventId", &m_eventId, buf_size, split_lvl);
  m_tree->Branch("size", &m_size, buf_size, split_lvl);
  m_tree->Branch("ipMomentumEst", &m_ipMomentumEst, buf_size, split_lvl);
  m_tree->Branch("vertexEst", &m_vertexEst, buf_size, split_lvl);
  m_tree->Branch("trackId", &m_trackId, buf_size, split_lvl);

  //------------------------------------------------------------------
  // Initialize the data handles
  m_seeds.initialize(m_cfg.inputSeeds);
}

ProcessCode RootSeedWriter::finalize() {
  if (m_file) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootSeedWriter::write(const AlgorithmContext& ctx) {
  const auto& inputSeeds = m_seeds(ctx);

  if (inputSeeds.empty()) {
    ACTS_DEBUG("Received empty seed vector. Continuing");
    return ProcessCode::SUCCESS;
  }

  std::lock_guard<std::mutex> lock(m_mutex);

  m_eventId = ctx.eventNumber;

  for (const auto& seed : inputSeeds) {
    const auto& sourceLinks =
        seed.sourceLinks | std::views::transform([](const auto& sl) {
          return sl.template get<SimpleSourceLink>();
        });
    if (sourceLinks.empty()) {
      continue;
    }

    m_seedMeasurementsGlob.clear();
    m_seedMeasurementsLoc.clear();
    m_geoIds.clear();

    m_seedMeasurementsGlob.reserve(sourceLinks.size());
    m_seedMeasurementsLoc.reserve(sourceLinks.size());
    m_geoIds.reserve(sourceLinks.size());
    for (const auto& sl : sourceLinks) {
      m_seedMeasurementsGlob.emplace_back(sl.parametersGlob().x(),
                                          sl.parametersGlob().y(),
                                          sl.parametersGlob().z());
      m_seedMeasurementsLoc.emplace_back(sl.parametersLoc().x(),
                                         sl.parametersLoc().y());
      m_geoIds.push_back(sl.geometryId().sensitive());
    }
    m_ipMomentumEst.SetPxPyPzE(
        seed.ipParameters.momentum().x(), seed.ipParameters.momentum().y(),
        seed.ipParameters.momentum().z(), seed.ipParameters.absoluteMomentum());
    m_vertexEst.SetXYZ(seed.ipParameters.position().x(),
                       seed.ipParameters.position().y(),
                       seed.ipParameters.position().z());

    m_size = sourceLinks.size();

    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}

namespace {

struct RootSeedWriterRegistrar {
  RootSeedWriterRegistrar() {
    using namespace TrackingPipeline;

    WriterRegistry::instance().registerBuilder(
      "RootSeedWriter",
      [](const toml::value& section,
         Acts::Logging::Level logLevel,
         const std::string& runRoot) -> WriterPtr {

        RootSeedWriter::Config cfg;

        cfg.inputSeeds =
            toml::find<std::string>(section, "inputSeeds");
        cfg.inputTruthClusters =
            toml::find<std::string>(section, "inputTruthClusters");
        cfg.treeName =
            toml::find<std::string>(section, "treeName");
        cfg.filePath =
            runRoot + "/" + toml::find<std::string>(section, "filePath");

        return std::make_shared<RootSeedWriter>(cfg, logLevel);
      });
  }
} _RootSeedWriterRegistrar;

}  // namespace
