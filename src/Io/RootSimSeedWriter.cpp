#include "TrackingPipeline/Io/RootSimSeedWriter.hpp"

#include "Acts/Utilities/Logger.hpp"
#include <Acts/Geometry/GeometryIdentifier.hpp>

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
  m_tree->Branch("size", &m_size, buf_size, split_lvl);
  m_tree->Branch("ipMomentumEst", &m_ipMomentumEst, buf_size, split_lvl);
  m_tree->Branch("vertexEst", &m_vertexEst, buf_size, split_lvl);

  // Truth info
  m_tree->Branch("trackId", &m_trackId, buf_size, split_lvl);
  m_tree->Branch("parentTrackId", &m_parentTrackId, buf_size, split_lvl);
  m_tree->Branch("runId", &m_runId, buf_size, split_lvl);
  m_tree->Branch("trueTrackSize", &m_trueTrackSize, buf_size, split_lvl);
  m_tree->Branch("trackInSeedSize", &m_trackInSeedSize, buf_size, split_lvl);
  m_tree->Branch("isSignal", &m_isSignal, buf_size, split_lvl);
  m_tree->Branch("ipMomentumTruth", &m_ipMomentumTruth, buf_size, split_lvl);
  m_tree->Branch("vertexTruth", &m_vertexTruth, buf_size, split_lvl);

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

  std::map<TrackID, std::set<Acts::GeometryIdentifier>> trackIds;
  for (const auto& cluster : inputTruthClusters) {
    if (!cluster.isSignal) {
      continue;
    }

    for (const auto& hit : cluster.truthHits) {
      trackIds[{hit.trackId, hit.parentTrackId, hit.runId}].insert(
          cluster.sourceLink.geometryId());
    }
  }

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

    std::map<TrackID, std::pair<int, std::set<Acts::GeometryIdentifier>>>
        seedSignalTrackIds;
    for (const auto& sl : sourceLinks) {
      const auto& cl = inputTruthClusters.at(sl.index());
      if (!cl.isSignal) {
        continue;
      }
      for (const auto& hit : cl.truthHits) {
        TrackID tid = {hit.trackId, hit.parentTrackId, hit.runId};
        seedSignalTrackIds[tid].first = sl.index();
        seedSignalTrackIds[tid].second.insert(sl.geometryId());
      }
    }
    if (seedSignalTrackIds.empty()) {
      m_trueTrackSize = 0;
      m_isSignal = false;
      m_ipMomentumTruth = TLorentzVector(0, 0, 0, 0);
      m_tree->Fill();
      continue;
    }

    auto maxTrack = std::max_element(
        seedSignalTrackIds.begin(), seedSignalTrackIds.end(),
        [](const auto& idA, const auto& idB) {
          return (idA.second.second.size() < idB.second.second.size());
        });
    const auto& [sigId, ic] = *maxTrack;
    auto [idx, geoIds] = ic;
    std::size_t trackSize = trackIds.at(sigId).size();

    std::tie(m_trackId, m_parentTrackId, m_runId) = sigId;

    m_trueTrackSize = trackSize;
    m_trackInSeedSize = geoIds.size();

    m_isSignal = true;

    const auto& cluster = inputTruthClusters.at(idx);
    for (const auto& hit : cluster.truthHits) {
      if (std::tie(hit.trackId, hit.parentTrackId, hit.runId) == sigId) {
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

  // Return success flag
  return ProcessCode::SUCCESS;
}
