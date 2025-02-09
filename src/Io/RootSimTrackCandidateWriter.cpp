#include "TrackingPipeline/Io/RootSimTrackCandidateWriter.hpp"

#include <algorithm>
#include <ranges>
#include <stdexcept>
#include <vector>

#include <bits/ranges_algo.h>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

RootSimTrackCandidateWriter::RootSimTrackCandidateWriter(
    const Config& config, Acts::Logging::Level level)
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
  int buf_size = 32000;
  int split_lvl = 0;

  m_tree->Branch("matchingDegree", &m_matchingDegree);
  m_tree->Branch("size", &m_size);
  m_tree->Branch("ipMomentumTruth", &m_ipMomentumTruth);

  //------------------------------------------------------------------
  // Initialize the data handles
  m_trackCandidates.initialize(m_cfg.inputTrackCandidates);
  m_truthClusters.initialize(m_cfg.inputTruthClusters);
}

ProcessCode RootSimTrackCandidateWriter::finalize() {
  if (m_file) {
    m_file->Write();
    m_file->Close();
    delete m_file;
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootSimTrackCandidateWriter::write(const AlgorithmContext& ctx) {
  auto inputTrackCandidates = m_trackCandidates(ctx);

  auto inputTruthClusters = m_truthClusters(ctx);

  std::lock_guard<std::mutex> lock(m_mutex);

  // Collect true track statistics
  std::map<TrackID, std::int32_t> trueTracksSig;

  // Collect true track statistics
  auto trueTrackIds =
      inputTruthClusters |
      std::views::filter([](const auto& cl) { return cl.isSignal; }) |
      std::views::transform([](const auto& cl) { return cl.truthHits; }) |
      std::views::join | std::views::transform([](const auto& hit) -> TrackID {
        return {hit.trackId, hit.parentTrackId, hit.runId};
      });

  for (const auto& id : trueTrackIds) {
    if (!trueTracksSig.contains(id)) {
      trueTracksSig[id] = 1;
    } else {
      trueTracksSig.at(id)++;
    }
  }

  for (const auto& candidate : inputTrackCandidates) {
    const auto& sourceLinks =
        candidate.sourceLinks | std::views::transform([](const auto& sl) {
          return sl.template get<SimpleSourceLink>();
        });
    m_size = sourceLinks.size();

    // Collect the truth info
    SimClusters candidateClusters;
    std::map<TrackID, std::vector<std::int32_t>> candidateTrackIds;
    candidateClusters.reserve(sourceLinks.size());
    for (const auto& sl : sourceLinks) {
      SimCluster cluster = inputTruthClusters.at(sl.index());
      candidateClusters.push_back(cluster);

      // Get the true hit
      TrackID currentTrackId;
      if (cluster.truthHits.size() == 0 || !cluster.isSignal) {
        currentTrackId = std::make_tuple(-1, -1, -1);
      } else {
        auto sig = std::ranges::find_if(cluster.truthHits, [](const auto& hit) {
          return (hit.trackId == 1);
        });

        currentTrackId =
            std::make_tuple(sig->trackId, sig->parentTrackId, sig->runId);
      }
      if (!candidateTrackIds.contains(currentTrackId)) {
        candidateTrackIds[currentTrackId] = {sl.index()};
      } else {
        candidateTrackIds.at(currentTrackId).push_back(sl.index());
      }
    }
    m_matchingDegree = 0;

    // Matching degree is computed with respect
    // to the most often occuring signal track
    auto refTrackId = std::ranges::max_element(
        candidateTrackIds, [](const auto& pairA, const auto& pairB) {
          if (std::get<0>(pairA.first) == -1 &&
              std::get<0>(pairB.first) != -1) {
            return true;
          } else if (std::get<0>(pairA.first) != -1 &&
                     std::get<0>(pairB.first) == -1) {
            return false;
          } else {
            return pairA.second.size() < pairB.second.size();
          }
        });
    if (std::get<0>(refTrackId->first) == -1) {
      m_matchingDegree = 0;
    } else {
      // Compute matching degree
      double trueTrackSize =
          std::ranges::count(trueTrackIds, refTrackId->first);

      if (trueTrackSize != m_cfg.targetTrueTrackSize) {
        continue;
      }

      m_matchingDegree = refTrackId->second.size() / trueTrackSize;
    }

    const auto& refCluster = inputTruthClusters.at(refTrackId->second.at(0));

    const auto& refHit = std::ranges::find_if(
        refCluster.truthHits, [&refTrackId](const auto& hit) {
          TrackID id = {hit.trackId, hit.parentTrackId, hit.runId};
          return id == refTrackId->first;
        });

    double me = 0.511 * Acts::UnitConstants::MeV;
    m_ipMomentumTruth.SetPxPyPzE(
        refHit->ipParameters.momentum().x(),
        refHit->ipParameters.momentum().y(),
        refHit->ipParameters.momentum().z(),
        std::hypot(refHit->ipParameters.absoluteMomentum(), me));

    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
