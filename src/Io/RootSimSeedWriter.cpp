#include "TrackingPipeline/Io/RootSimSeedWriter.hpp"

#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>

#include <cstddef>
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
  int bufSize = 32000;
  int splitLvl = 0;

  // Measurements
  m_tree->Branch("measurementsGlob", &m_seedMeasurementsGlob, bufSize,
                 splitLvl);
  m_tree->Branch("measurementsLoc", &m_seedMeasurementsLoc, bufSize, splitLvl);
  m_tree->Branch("geoIds", &m_geoIds, bufSize, splitLvl);

  // Seed properties
  m_tree->Branch("eventId", &m_eventId, bufSize, splitLvl);
  m_tree->Branch("size", &m_size, bufSize, splitLvl);
  m_tree->Branch("originMomentumEst", &m_originMomentumEst, bufSize, splitLvl);
  m_tree->Branch("vertexEst", &m_vertexEst, bufSize, splitLvl);

  // Truth info
  m_tree->Branch("trackId", &m_trackId, bufSize, splitLvl);
  m_tree->Branch("parentTrackId", &m_parentTrackId, bufSize, splitLvl);
  m_tree->Branch("runId", &m_runId, bufSize, splitLvl);
  m_tree->Branch("trueTrackSize", &m_trueTrackSize, bufSize, splitLvl);
  m_tree->Branch("trackInSeedSize", &m_trackInSeedSize, bufSize, splitLvl);
  m_tree->Branch("isSignal", &m_isSignal, bufSize, splitLvl);
  m_tree->Branch("originMomentumTruth", &m_originMomentumTruth, bufSize,
                 splitLvl);
  m_tree->Branch("vertexTruth", &m_vertexTruth, bufSize, splitLvl);

  //------------------------------------------------------------------
  // Initialize the data handles
  m_inputSeeds.initialize(m_cfg.inputSeeds);
  m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
  m_inputSimClusters.initialize(m_cfg.inputSimClusters);
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
}

ProcessCode RootSimSeedWriter::finalize() {
  if (m_file != nullptr) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootSimSeedWriter::write(const AlgorithmContext& ctx) {
  const auto& inputSeeds = m_inputSeeds(ctx);
  const auto& inputTrackParameters = m_inputTrackParameters(ctx);
  const auto& inputSimClusters = m_inputSimClusters(ctx);
  const auto& inputSourceLinks = m_inputSourceLinks(ctx);

  if (inputSeeds.empty()) {
    ACTS_DEBUG("Received empty seed vector. Skipping");
    return ProcessCode::SUCCESS;
  }

  std::lock_guard<std::mutex> lock(m_mutex);

  m_eventId = ctx.eventNumber;

  std::map<TrackID, std::set<Acts::GeometryIdentifier>> trackIds;
  for (const auto& cluster : inputSimClusters) {
    if (!cluster.isSignal) {
      continue;
    }

    for (const auto& hit : cluster.truthHits) {
      trackIds[{hit.trackId, hit.parentTrackId, hit.runId}].insert(
          cluster.sourceLink.get<SimpleSourceLink>().geometryId());
    }
  }

  for (const auto& seed : inputSeeds) {
    const auto& sourceLinkIndices = seed.sourceLinkIndices;
    m_size = sourceLinkIndices.size();

    m_seedMeasurementsGlob.clear();
    m_seedMeasurementsLoc.clear();
    m_geoIds.clear();

    m_seedMeasurementsGlob.reserve(m_size);
    m_seedMeasurementsLoc.reserve(m_size);
    m_geoIds.reserve(m_size);
    for (std::size_t idx : sourceLinkIndices) {
      const auto& ssl = inputSourceLinks.at(idx).get<SimpleSourceLink>();

      const Acts::Vector3& parsGlob = ssl.parametersGlob();
      m_seedMeasurementsGlob.emplace_back(parsGlob.x(), parsGlob.y(),
                                          parsGlob.z());

      const Acts::Vector2& parsLoc = ssl.parametersLoc();
      m_seedMeasurementsLoc.emplace_back(parsLoc.x(), parsLoc.y());
      m_geoIds.push_back(ssl.geometryId().sensitive());
    }

    const auto& originParameters =
        inputTrackParameters.at(seed.originParametersIndex);
    m_originMomentumEst.SetPxPyPzE(
        originParameters.momentum().x(), originParameters.momentum().y(),
        originParameters.momentum().z(), originParameters.absoluteMomentum());
    m_vertexEst.SetXYZ(originParameters.position().x(),
                       originParameters.position().y(),
                       originParameters.position().z());

    std::map<TrackID,
             std::pair<int, std::unordered_set<Acts::GeometryIdentifier>>>
        seedSignalTrackIds;
    for (std::size_t idx : sourceLinkIndices) {
      const auto& ssl = inputSourceLinks.at(idx).get<SimpleSourceLink>();
      const auto& cl = inputSimClusters.at(ssl.index());
      if (!cl.isSignal) {
        continue;
      }
      for (const auto& hit : cl.truthHits) {
        TrackID tid = {hit.trackId, hit.parentTrackId, hit.runId};
        seedSignalTrackIds[tid].first = ssl.index();
        seedSignalTrackIds[tid].second.insert(ssl.geometryId());
      }
    }
    if (seedSignalTrackIds.empty()) {
      m_trueTrackSize = 0;
      m_isSignal = false;
      m_originMomentumTruth = TLorentzVector(0, 0, 0, 0);
      m_tree->Fill();
      continue;
    }

    auto maxTrack = std::max_element(
        seedSignalTrackIds.begin(), seedSignalTrackIds.end(),
        [](const auto& idA, const auto& idB) {
          return (idA.second.second.size() < idB.second.second.size());
        });
    const auto& [sigId, ic] = *maxTrack;
    auto [sourceLinkIdx, geoIds] = ic;
    std::size_t trackSize = trackIds.at(sigId).size();

    std::tie(m_trackId, m_parentTrackId, m_runId) = sigId;

    m_trueTrackSize = trackSize;
    m_trackInSeedSize = geoIds.size();

    m_isSignal = true;

    const auto& cluster = inputSimClusters.at(sourceLinkIdx);
    for (const auto& hit : cluster.truthHits) {
      if (std::tie(hit.trackId, hit.parentTrackId, hit.runId) == sigId) {
        const Acts::Vector3& momentum = hit.ipParameters.momentum();
        const Acts::Vector3& position = hit.ipParameters.position();
        double absMom = hit.ipParameters.absoluteMomentum();

        m_originMomentumTruth.SetPxPyPzE(momentum.x(), momentum.y(),
                                         momentum.z(), absMom);
        m_vertexTruth.SetXYZ(position.x(), position.y(), position.z());
        break;
      }
    }

    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
