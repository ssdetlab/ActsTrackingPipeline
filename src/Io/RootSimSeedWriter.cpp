#include "TrackingPipeline/Io/RootSimSeedWriter.hpp"

#include <algorithm>
#include <ranges>
#include <vector>

#include <TLorentzVector.h>
#include <bits/ranges_algo.h>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

RootSimSeedWriter::RootSimSeedWriter(const Config& config,
                                     Acts::Logging::Level level)
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
  m_tree->Branch("efficiency", &m_efficiency);
  m_tree->Branch("size", &m_size);
  m_tree->Branch("isSignal", &m_isSignal);
  m_tree->Branch("ipMomentumTruth", &m_ipMomentumTruth);

  //------------------------------------------------------------------
  // Initialize the data handles
  m_seeds.initialize(m_cfg.inputSeeds);
  m_truthClusters.initialize(m_cfg.inputTruthClusters);
}

ProcessCode RootSimSeedWriter::finalize() {
  if (m_file) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootSimSeedWriter::write(const AlgorithmContext& ctx) {
  auto inputSeeds = m_seeds(ctx);

  auto inputTruthClusters = m_truthClusters(ctx);

  std::lock_guard<std::mutex> lock(m_mutex);

  // Collect true track statistics
  auto trackIds =
      inputTruthClusters |
      std::views::filter([](const auto& cl) { return cl.isSignal; }) |
      std::views::transform([](const auto& cl) { return cl.truthHits; }) |
      std::views::join | std::views::transform([](const auto& hit) -> TrackID {
        return {hit.trackId, hit.parentTrackId, hit.runId};
      });

  for (const auto& seed : inputSeeds) {
    const auto& sourceLinks =
        seed.sourceLinks | std::views::transform([](const auto& sl) {
          return sl.template get<SimpleSourceLink>();
        });
    m_size = sourceLinks.size();
    if (!inputTruthClusters.at(sourceLinks.front().index()).isSignal) {
      m_efficiency = 0;
      m_isSignal = false;
      m_ipMomentumTruth = TLorentzVector(0, 0, 0, 0);
      m_tree->Fill();
      continue;
    }

    SimClusters seedClusters;
    seedClusters.reserve(sourceLinks.size());
    for (const auto& sl : sourceLinks) {
      seedClusters.push_back(inputTruthClusters.at(sl.index()));
    }

    auto pivotCluster = seedClusters.front();

    auto signalTrackIds =
        pivotCluster.truthHits |
        std::views::filter([](const auto& hit) { return (hit.trackId == 1); }) |
        std::views::transform([](const auto& hit) -> TrackID {
          return {hit.trackId, hit.parentTrackId, hit.runId};
        });

    auto seedTrackIds =
        seedClusters |
        std::views::filter([](const auto& cl) { return cl.isSignal; }) |
        std::views::transform([](const auto& cl) { return cl.truthHits; }) |
        std::views::join |
        std::views::transform([](const auto& hit) -> TrackID {
          return {hit.trackId, hit.parentTrackId, hit.runId};
        });

    for (const auto& sigId : signalTrackIds) {
      double trackSize = std::ranges::count(trackIds, sigId);
      double trackInSeedSize = std::ranges::count(seedTrackIds, sigId);

      if (trackSize != m_cfg.targetTrueTrackSize) {
        continue;
      }

      m_efficiency = trackInSeedSize / trackSize;

      m_isSignal = true;
      double me = 0.511 * Acts::UnitConstants::MeV;
      for (const auto& hit : pivotCluster.truthHits) {
        if (hit.trackId == 1) {
          m_ipMomentumTruth.SetPxPyPzE(
              hit.ipParameters.momentum().x(), hit.ipParameters.momentum().y(),
              hit.ipParameters.momentum().z(),
              std::hypot(hit.ipParameters.absoluteMomentum(), me));

          break;
        }
      }

      m_tree->Fill();
    }
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
