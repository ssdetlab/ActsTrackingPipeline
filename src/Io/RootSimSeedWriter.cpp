#include "TrackingPipeline/Io/RootSimSeedWriter.hpp"
#include <Acts/Utilities/Logger.hpp>

#include <ranges>
#include <vector>

#include "TLorentzVector.h"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

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
  m_tree->Branch("size", &m_size);
  m_tree->Branch("ipMomentumEst", &m_ipMomentumEst);
  m_tree->Branch("vertexEst", &m_vertexEst);

  // Truth info
  m_tree->Branch("trueTrackSize", &m_trueTrackSize);
  m_tree->Branch("trackInSeedSize", &m_trackInSeedSize);
  m_tree->Branch("matchingDegree", &m_matchingDegree);
  m_tree->Branch("isSignal", &m_isSignal);
  m_tree->Branch("ipMomentumTruth", &m_ipMomentumTruth);
  m_tree->Branch("vertexTruth", &m_vertexTruth);

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
  const auto& inputSeeds = m_seeds(ctx);
  const auto& inputTruthClusters = m_truthClusters(ctx);

  if (inputSeeds.empty()) {
    ACTS_DEBUG("Received empty seed vector. Continuing");
    return ProcessCode::SUCCESS;
  }

  std::lock_guard<std::mutex> lock(m_mutex);

  m_eventId = ctx.eventNumber;

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

    if (!inputTruthClusters.at(sourceLinks.front().index()).isSignal) {
      m_trueTrackSize = 0;
      m_matchingDegree = 0;
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

    const auto& pivotCluster = seedClusters.front();

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
      std::size_t trackSize = std::ranges::count(trackIds, sigId);
      std::size_t trackInSeedSize = std::ranges::count(seedTrackIds, sigId);

      m_trueTrackSize = trackSize;
      m_trackInSeedSize = trackInSeedSize;
      m_matchingDegree = trackInSeedSize / static_cast<double>(trackSize);

      m_isSignal = true;
      // Note: This assumes that there's only 1 track in a cluster
      for (const auto& hit : pivotCluster.truthHits) {
        if (hit.trackId == 1) {
          m_ipMomentumTruth.SetPxPyPzE(hit.ipParameters.momentum().x(),
                                       hit.ipParameters.momentum().y(),
                                       hit.ipParameters.momentum().z(),
                                       hit.ipParameters.absoluteMomentum());
          m_vertexTruth.SetXYZ(hit.ipParameters.position().x(),
                               hit.ipParameters.position().y(),
                               hit.ipParameters.position().z());

          break;
        }
      }

      m_tree->Fill();
    }
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
