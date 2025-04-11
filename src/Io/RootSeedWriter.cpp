#include "TrackingPipeline/Io/RootSeedWriter.hpp"

#include <ranges>
#include <vector>

#include <TLorentzVector.h>
#include <bits/ranges_algo.h>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

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
  // Track tree branches
  m_tree->Branch("size", &m_size);

  //------------------------------------------------------------------
  // Initialize the data handles
  m_seeds.initialize(m_cfg.inputSeeds);
}

ProcessCode RootSeedWriter::finalize() {
  if (m_file) {
    m_file->Write();
    m_file->Close();
    delete m_file;
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootSeedWriter::write(const AlgorithmContext& ctx) {
  auto inputSeeds = m_seeds(ctx);

  std::lock_guard<std::mutex> lock(m_mutex);

  for (const auto& seed : inputSeeds) {
    const auto& sourceLinks =
        seed.sourceLinks | std::views::transform([](const auto& sl) {
          return sl.template get<SimpleSourceLink>();
        });
    m_size = sourceLinks.size();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
