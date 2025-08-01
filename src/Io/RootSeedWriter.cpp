#include "TrackingPipeline/Io/RootSeedWriter.hpp"

#include <ranges>
#include <vector>

#include <bits/ranges_algo.h>

#include "TLorentzVector.h"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

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
  m_tree->Branch("pivot", &m_geoCenterPivot, buf_size, split_lvl);
  m_tree->Branch("measurements", &m_seedMeasurements, buf_size, split_lvl);

  // Seed properties
  m_tree->Branch("nLayes", &m_nLayers);
  m_tree->Branch("size", &m_size);
  m_tree->Branch("ipMomentumEst", &m_ipMomentumEst);
  m_tree->Branch("vertexEst", &m_vertexEst);

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
  auto inputSeeds = m_seeds(ctx);

  std::lock_guard<std::mutex> lock(m_mutex);

  double me = 0.511 * Acts::UnitConstants::MeV;
  for (const auto& seed : inputSeeds) {
    std::set<Acts::GeometryIdentifier> seedGeoIds;

    const auto& sourceLinks =
        seed.sourceLinks | std::views::transform([](const auto& sl) {
          return sl.template get<SimpleSourceLink>();
        });

    const auto& pivotSl = sourceLinks.front();
    const auto* surfPivot = m_cfg.surfaceAccessor(Acts::SourceLink(pivotSl));
    Acts::Vector3 geoCenterGlobal = surfPivot->localToGlobal(
        ctx.geoContext, pivotSl.parameters(), Acts::Vector3(0, 1, 0));
    m_geoCenterPivot.SetX(geoCenterGlobal.x());
    m_geoCenterPivot.SetY(geoCenterGlobal.y());
    m_geoCenterPivot.SetZ(geoCenterGlobal.z());

    m_seedMeasurements.clear();
    m_seedMeasurements.reserve(sourceLinks.size());
    for (const auto& sl : sourceLinks) {
      const auto* surf = m_cfg.surfaceAccessor(Acts::SourceLink(sl));
      seedGeoIds.insert(surf->geometryId());
      Acts::Vector3 measGlobal = surf->localToGlobal(
          ctx.geoContext, sl.parameters(), Acts::Vector3(0, 1, 0));
      m_seedMeasurements.emplace_back(measGlobal.x(), measGlobal.y(),
                                      measGlobal.z());
    }
    m_ipMomentumEst.SetPxPyPzE(
        seed.ipParameters.momentum().x(), seed.ipParameters.momentum().y(),
        seed.ipParameters.momentum().z(),
        std::hypot(seed.ipParameters.absoluteMomentum(), me));
    m_size = sourceLinks.size();
    m_nLayers = seedGeoIds.size();

    m_vertexEst.SetX(seed.ipParameters.position().x());
    m_vertexEst.SetY(seed.ipParameters.position().y());
    m_vertexEst.SetZ(seed.ipParameters.position().z());

    m_tree->Fill();
  }
  return ProcessCode::SUCCESS;
}
