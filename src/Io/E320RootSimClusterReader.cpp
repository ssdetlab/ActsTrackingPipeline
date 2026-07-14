#include "TrackingPipeline/Io/E320RootSimClusterReader.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include <Acts/Definitions/PdgParticle.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/ExtendedSourceLink.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

using namespace Acts::UnitLiterals;

namespace E320 {

E320RootSimClusterReader::E320RootSimClusterReader(const Config& config,
                                                   Acts::Logging::Level level)
    : IReader(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.filePaths.empty()) {
    throw std::invalid_argument("Missing filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  if (m_cfg.filePaths.size() == 1) {
    m_file = new TFile(m_cfg.filePaths.at(0).c_str());
    m_tree = m_file->Get<TTree>(m_cfg.treeName.c_str());
  } else {
    m_chainOwner = new TChain(m_cfg.treeName.c_str());
    // Add the files to the chain
    for (const auto& path : m_cfg.filePaths) {
      m_chainOwner->Add(path.c_str());
    }
    m_tree = dynamic_cast<TTree*>(m_chainOwner);
  }

  //------------------------------------------------------------------
  // Tree branches
  int bufSize = 32000;
  int splitLvl = 0;

  // Cluster parameters
  m_tree->SetBranchAddress("geoCenterGlobal", &m_geoCenterGlobal);
  m_tree->SetBranchAddress("geoCenterLocal", &m_geoCenterLocal);
  m_tree->SetBranchAddress("clusterCov", &m_clusterCov);
  m_tree->SetBranchAddress("geoId", &m_geoId);

  // Measurement hits
  m_tree->SetBranchAddress("trackHitsGlobal", &m_trackHitsGlobal);
  m_tree->SetBranchAddress("trackHitsLocal", &m_trackHitsLocal);
  m_tree->SetBranchAddress("eventId", &m_eventId);
  m_tree->SetBranchAddress("charge", &m_charge);
  m_tree->SetBranchAddress("pdgId", &m_pdgId);

  // Bound origin parameters
  m_tree->SetBranchAddress("boundTrackParameters", &m_boundTrackParameters);
  m_tree->SetBranchAddress("boundTrackCov", &m_boundTrackCov);

  // Origin momentum
  m_tree->SetBranchAddress("originMomentum", &m_originMomentum);

  // Origin vertex
  m_tree->SetBranchAddress("vertex", &m_vertex);

  // Momentum at clusters
  m_tree->SetBranchAddress("onSurfaceMomentum", &m_onSurfaceMomentum);

  // Track ID
  m_tree->SetBranchAddress("trackId", &m_trackId);
  m_tree->SetBranchAddress("parentTrackId", &m_parentTrackId);
  m_tree->SetBranchAddress("runId", &m_runId);

  // Misc
  m_tree->SetBranchAddress("isSignal", &m_isSignal);

  // Disable all branches and only enable event-id for a first scan of the
  // file
  m_tree->SetBranchStatus("*", false);
  if (m_tree->GetBranch("eventId") == nullptr) {
    throw std::invalid_argument("Missing eventId branch");
  }
  m_tree->SetBranchStatus("eventId", true);
  auto nEntries = static_cast<std::size_t>(m_tree->GetEntries());

  // Go through all entries and store the position of the events
  m_tree->GetEntry(0);
  m_eventMap.emplace_back(m_eventId, 0, 0);
  for (std::size_t i = 0; i < nEntries; ++i) {
    m_tree->GetEntry(i);
    if (m_eventId != std::get<0>(m_eventMap.back())) {
      std::get<2>(m_eventMap.back()) = i;
      m_eventMap.emplace_back(m_eventId, i, i);
    }
    if (i == nEntries - 1) {
      std::get<2>(m_eventMap.back()) = nEntries;
    }
  }

  // Sort by event id
  std::ranges::sort(m_eventMap, [](const auto& a, const auto& b) {
    return std::get<0>(a) < std::get<0>(b);
  });

  // Re-Enable all branches
  m_tree->SetBranchStatus("*", true);
  ACTS_DEBUG("Event range: " << availableEvents().first << " - "
                             << availableEvents().second);

  // Initialize BPM geometry ids
  const auto& goInst = *GeometryOptions::instance();
  m_bpmGeoIds.insert(goInst.bpm0Parameters.geoId);
  m_bpmGeoIds.insert(goInst.bpm1Parameters.geoId);
  m_bpmGeoIds.insert(goInst.bpm2Parameters.geoId);
  m_bpmGeoIds.insert(goInst.bpm3Parameters.geoId);
  m_bpmGeoIds.insert(goInst.ipSurfaceParameters.geoId);

  // Initialize the data handles
  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
  m_outputSimClusters.initialize(m_cfg.outputSimClusters);
  m_outputDetSourceLinkIndices.initialize(m_cfg.outputDetSourceLinkIndices);
  m_outputBpmSourceLinkIndices.initialize(m_cfg.outputBpmSourceLinkIndices);
}

std::pair<std::size_t, std::size_t> E320RootSimClusterReader::availableEvents()
    const {
  return {std::get<0>(m_eventMap.front()), std::get<0>(m_eventMap.back()) + 1};
}

ProcessCode E320RootSimClusterReader::read(const AlgorithmContext& ctx) {
  auto it = std::ranges::find_if(m_eventMap, [&](const auto& a) {
    return std::get<0>(a) == ctx.eventNumber;
  });

  if (it == m_eventMap.end()) {
    // Explicitly warn if it happens for the first or last event as that might
    // indicate a human error
    if ((ctx.eventNumber == availableEvents().first) &&
        (ctx.eventNumber == availableEvents().second - 1)) {
      ACTS_WARNING("Reading empty event: " << ctx.eventNumber);
    } else {
      ACTS_DEBUG("Reading empty event: " << ctx.eventNumber);
    }

    m_outputSourceLinks(ctx, {});
    m_outputSimClusters(ctx, {});
    m_outputDetSourceLinkIndices(ctx, {});
    m_outputBpmSourceLinkIndices(ctx, {});

    // Return success flag
    return ProcessCode::SUCCESS;
  }

  // lock the mutex
  std::lock_guard<std::mutex> lock(m_read_mutex);

  ACTS_DEBUG("Reading event: " << std::get<0>(*it)
                               << " stored in entries: " << std::get<1>(*it)
                               << " - " << std::get<2>(*it));

  // Create the measurements
  std::vector<Acts::SourceLink> sourceLinks{};
  std::vector<std::size_t> detSourceLinkIndices{};
  std::vector<std::size_t> bpmSourceLinkIndices{};
  SimClusters simClusters{};
  std::size_t eventId = std::get<0>(*it);
  for (auto entry = std::get<1>(*it); entry < std::get<2>(*it); entry++) {
    m_tree->GetEntry(entry);

    if (m_geoId < m_cfg.minGeoId || m_geoId > m_cfg.maxGeoId) {
      continue;
    }

    Acts::GeometryIdentifier geoId;
    geoId.setSensitive(m_geoId);

    Acts::Vector2 geoCenterLocal(m_geoCenterLocal->X(), m_geoCenterLocal->Y());

    Acts::Vector3 geoCenterGlobal;
    if (m_cfg.surfaceLocalToGlobal) {
      geoCenterGlobal = m_cfg.surfaceMap.at(geoId)->localToGlobal(
          ctx.geoContext, geoCenterLocal, Acts::Vector3::UnitX());
    } else {
      geoCenterGlobal =
          Acts::Vector3(m_geoCenterGlobal->X(), m_geoCenterGlobal->Y(),
                        m_geoCenterGlobal->Z());
    }

    Acts::SquareMatrix2 clusterCov;
    clusterCov << (*m_clusterCov)(0, 0), (*m_clusterCov)(0, 1),
        (*m_clusterCov)(1, 0), (*m_clusterCov)(1, 1);

    SimpleSourceLink obsSourceLink(geoCenterLocal, geoCenterGlobal, clusterCov,
                                   geoId, eventId, sourceLinks.size());
    if (m_bpmGeoIds.contains(geoId.sensitive())) {
      bpmSourceLinkIndices.push_back(sourceLinks.size());
    } else {
      detSourceLinkIndices.push_back(sourceLinks.size());
    }
    sourceLinks.push_back(Acts::SourceLink(obsSourceLink));

    SimHits hits;
    hits.reserve(m_trackHitsLocal->size());
    for (std::size_t i = 0; i < m_trackHitsLocal->size(); i++) {
      Acts::Vector4 vertex(m_vertex->at(i).X(), m_vertex->at(i).Y(),
                           m_vertex->at(i).Z(), 0);
      Acts::Vector3 ipDirection(m_originMomentum->at(i).X(),
                                m_originMomentum->at(i).Y(),
                                m_originMomentum->at(i).Z());
      ipDirection.normalize();

      Acts::BoundMatrix ipCov = Acts::BoundMatrix::Zero();
      for (std::size_t n = 0; n < Acts::eBoundSize; n++) {
        for (std::size_t m = 0; m < Acts::eBoundSize; m++) {
          ipCov(n, m) = m_boundTrackCov->at(i)(n, m);
        }
      }

      Acts::ParticleHypothesis hypothesis(Acts::PdgParticle(m_pdgId->at(i)));
      Acts::CurvilinearTrackParameters ipParameters(
          vertex, ipDirection, m_charge->at(i) / m_originMomentum->at(i).P(),
          ipCov, hypothesis);

      Acts::BoundVector truthParameters;
      truthParameters[Acts::eBoundLoc0] = m_trackHitsLocal->at(i).X();
      truthParameters[Acts::eBoundLoc1] = m_trackHitsLocal->at(i).Y();
      truthParameters[Acts::eBoundPhi] = m_onSurfaceMomentum->at(i).Phi();
      truthParameters[Acts::eBoundTheta] = m_onSurfaceMomentum->at(i).Theta();
      truthParameters[Acts::eBoundQOverP] =
          m_charge->at(i) / m_onSurfaceMomentum->at(i).P();

      Acts::Vector3 trueTrackHitGlobal(m_trackHitsGlobal->at(i).X(),
                                       m_trackHitsGlobal->at(i).Y(),
                                       m_trackHitsGlobal->at(i).Z());

      SimHit trackHit{std::move(truthParameters), std::move(trueTrackHitGlobal),
                      std::move(ipParameters),    m_trackId->at(i),
                      m_parentTrackId->at(i),     m_runId->at(i)};
      hits.push_back(trackHit);
    }
    SimCluster cluster{obsSourceLink, hits, static_cast<bool>(m_isSignal)};
    simClusters.push_back(cluster);
  }

  ACTS_DEBUG("Read " << sourceLinks.size() << " source links");
  ACTS_DEBUG("Read " << simClusters.size() << " clusters");
  ACTS_DEBUG("Read " << detSourceLinkIndices.size()
                     << " detector source link indices");
  ACTS_DEBUG("Read " << bpmSourceLinkIndices.size()
                     << " BPM source link indices");

  m_outputSourceLinks(ctx, std::move(sourceLinks));
  m_outputSimClusters(ctx, std::move(simClusters));
  m_outputDetSourceLinkIndices(ctx, std::move(detSourceLinkIndices));
  m_outputBpmSourceLinkIndices(ctx, std::move(bpmSourceLinkIndices));

  // Return success flag
  return ProcessCode::SUCCESS;
}

}  // namespace E320
