#include "TrackingPipeline/Io/RootSimClusterReader.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include <Acts/Definitions/PdgParticle.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <stdexcept>
#include <vector>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

using namespace Acts::UnitLiterals;

RootSimClusterReader::RootSimClusterReader(const Config& config,
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

  m_chain = new TChain(m_cfg.treeName.c_str());

  //------------------------------------------------------------------
  // Tree branches
  int buf_size = 32000;
  int split_lvl = 0;

  // Parameters at measurements
  m_chain->SetBranchAddress("geoCenterLocal", &m_geoCenterLocal);
  m_chain->SetBranchAddress("geoCenterGlobal", &m_geoCenterGlobal);
  m_chain->SetBranchAddress("cov", &m_cov);
  m_chain->SetBranchAddress("geoId", &m_geoId);
  m_chain->SetBranchAddress("trackHitsLocal", &m_trackHitsLoc);
  m_chain->SetBranchAddress("trackHitsGlobal", &m_trackHitsGlob);
  m_chain->SetBranchAddress("eventId", &m_eventId);
  m_chain->SetBranchAddress("charge", &m_charge);
  m_chain->SetBranchAddress("pdgId", &m_pdgId);

  // Parameters at the origin
  m_chain->SetBranchAddress("originMomentum", &m_originMomentum);
  m_chain->SetBranchAddress("vertex", &m_vertex);
  m_chain->SetBranchAddress("onSurfaceMomentum", &m_onSurfaceMomentum);

  // Track ID
  m_chain->SetBranchAddress("trackId", &m_trackId);
  m_chain->SetBranchAddress("parentTrackId", &m_parentTrackId);
  m_chain->SetBranchAddress("runId", &m_runId);

  // Misc
  m_chain->SetBranchAddress("isSignal", &m_isSignal);

  // Add the files to the chain
  for (const auto& path : m_cfg.filePaths) {
    m_chain->Add(path.c_str());
  }

  // Disable all branches and only enable event-id for a first scan of the
  // file
  m_chain->SetBranchStatus("*", false);
  if (!m_chain->GetBranch("eventId")) {
    throw std::invalid_argument("Missing eventId branch");
  }
  m_chain->SetBranchStatus("eventId", true);
  auto nEntries = static_cast<std::size_t>(m_chain->GetEntries());

  // Go through all entries and store the position of the events
  m_chain->GetEntry(0);
  m_eventMap.emplace_back(m_eventId, 0, 0);
  for (std::size_t i = 0; i < nEntries; ++i) {
    m_chain->GetEntry(i);
    if (m_eventId != std::get<0>(m_eventMap.back())) {
      std::get<2>(m_eventMap.back()) = i;
      m_eventMap.emplace_back(m_eventId, i, i);
    }
  }

  // Sort by event id
  std::ranges::sort(m_eventMap, [](const auto& a, const auto& b) {
    return std::get<0>(a) < std::get<0>(b);
  });

  std::get<2>(m_eventMap.back()) = nEntries;

  // Re-Enable all branches
  m_chain->SetBranchStatus("*", true);
  ACTS_DEBUG("Event range: " << availableEvents().first << " - "
                             << availableEvents().second);

  // Initialize the data handles
  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
  m_outputSimClusters.initialize(m_cfg.outputSimClusters);

  Acts::BoundVector ipStdDev;
  ipStdDev[Acts::eBoundLoc0] = 100_um;
  ipStdDev[Acts::eBoundLoc1] = 100_um;
  ipStdDev[Acts::eBoundTime] = 25_ns;
  ipStdDev[Acts::eBoundPhi] = 2_degree;
  ipStdDev[Acts::eBoundTheta] = 2_degree;
  ipStdDev[Acts::eBoundQOverP] = 1 / 100_GeV;
  m_ipCov = ipStdDev.cwiseProduct(ipStdDev).asDiagonal();
}

std::pair<std::size_t, std::size_t> RootSimClusterReader::availableEvents()
    const {
  return {std::get<0>(m_eventMap.front()), std::get<0>(m_eventMap.back()) + 1};
}

ProcessCode RootSimClusterReader::read(const AlgorithmContext& ctx) {
  auto it = std::ranges::find_if(m_eventMap, [&](const auto& a) {
    return std::get<0>(a) == ctx.eventNumber;
  });

  if (it == m_eventMap.end()) {
    // explicitly warn if it happens for the first or last event as that might
    // indicate a human error
    if ((ctx.eventNumber == availableEvents().first) &&
        (ctx.eventNumber == availableEvents().second - 1)) {
      ACTS_WARNING("Reading empty event: " << ctx.eventNumber);
    } else {
      ACTS_DEBUG("Reading empty event: " << ctx.eventNumber);
    }

    m_outputSourceLinks(ctx, {});
    m_outputSimClusters(ctx, {});

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
  SimClusters simClusters{};
  std::size_t eventId = std::get<0>(*it);
  std::size_t sslIdx = 0;
  for (auto entry = std::get<1>(*it); entry < std::get<2>(*it); entry++) {
    m_chain->GetEntry(entry);

    // if (m_geoId < 20) {
    //   continue;
    // }

    Acts::Vector2 geoCenterLocal(m_geoCenterLocal->X(), m_geoCenterLocal->Y());
    Acts::Vector3 geoCenterGlobal(
        m_geoCenterGlobal->X(), m_geoCenterGlobal->Y(), m_geoCenterGlobal->Z());
    Acts::ActsSquareMatrix<2> cov;
    cov << (*m_cov)(0, 0), (*m_cov)(0, 1), (*m_cov)(1, 0), (*m_cov)(1, 1);

    Acts::GeometryIdentifier geoId;
    geoId.setSensitive(m_geoId);
    SimpleSourceLink obsSourceLink(geoCenterLocal, geoCenterGlobal, cov, geoId,
                                   eventId, sslIdx);
    sourceLinks.push_back(Acts::SourceLink{obsSourceLink});

    SimHits hits;
    hits.reserve(m_trackHitsLoc->size());
    for (std::size_t i = 0; i < m_trackHitsLoc->size(); i++) {
      Acts::Vector4 vertex(m_vertex->at(i).X(), m_vertex->at(i).Y(),
                           m_vertex->at(i).Z(), 0);
      Acts::Vector3 ipDirection(m_originMomentum->at(i).X(),
                                m_originMomentum->at(i).Y(),
                                m_originMomentum->at(i).Z());
      ipDirection.normalize();
      Acts::ParticleHypothesis hypothesis(Acts::PdgParticle(m_pdgId->at(i)));
      Acts::CurvilinearTrackParameters ipParameters(
          vertex, ipDirection, m_charge->at(i) / m_originMomentum->at(i).P(),
          m_ipCov, hypothesis);

      Acts::BoundVector truthParameters;
      truthParameters[Acts::eBoundLoc0] = m_trackHitsLoc->at(i).X();
      truthParameters[Acts::eBoundLoc1] = m_trackHitsLoc->at(i).Y();
      truthParameters[Acts::eBoundPhi] = m_onSurfaceMomentum->at(i).Phi();
      truthParameters[Acts::eBoundTheta] = m_onSurfaceMomentum->at(i).Theta();
      truthParameters[Acts::eBoundQOverP] =
          m_charge->at(i) / m_onSurfaceMomentum->at(i).P();

      Acts::Vector3 trueTrackHitGlobal(m_trackHitsGlob->at(i).X(),
                                       m_trackHitsGlob->at(i).Y(),
                                       m_trackHitsGlob->at(i).Z());

      SimHit trackHit{std::move(truthParameters), std::move(trueTrackHitGlobal),
                      std::move(ipParameters),    m_trackId->at(i),
                      m_parentTrackId->at(i),     m_runId->at(i)};
      hits.push_back(trackHit);
    }
    SimCluster cluster{obsSourceLink, hits, static_cast<bool>(m_isSignal)};
    simClusters.push_back(cluster);
    sslIdx++;
  }

  ACTS_DEBUG("Read " << sourceLinks.size() << " source links");
  ACTS_DEBUG("Read " << simClusters.size() << " clusters");
  m_outputSourceLinks(ctx, std::move(sourceLinks));
  m_outputSimClusters(ctx, std::move(simClusters));

  // Return success flag
  return ProcessCode::SUCCESS;
}
