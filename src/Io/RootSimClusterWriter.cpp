#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

RootSimClusterWriter::RootSimClusterWriter(const Config& config,
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

  // Parameters at measurements
  m_tree->Branch("geoCenterGlobal", &m_geoCenterGlobal, buf_size, split_lvl);
  m_tree->Branch("geoCenterLocal", &m_geoCenterLocal, buf_size, split_lvl);
  m_tree->Branch("cov", &m_cov, buf_size, split_lvl);
  m_tree->Branch("geoId", &m_geoId, buf_size, split_lvl);
  m_tree->Branch("trackHitsGlobal", &m_trackHitsGlobal, buf_size, split_lvl);
  m_tree->Branch("trackHitsLocal", &m_trackHitsLocal, buf_size, split_lvl);
  m_tree->Branch("eventId", &m_eventId, buf_size, split_lvl);
  m_tree->Branch("charge", &m_charge, buf_size, split_lvl);
  m_tree->Branch("pdgId", &m_pdgId, buf_size, split_lvl);

  // Parameters at the origin
  m_tree->Branch("originMomentum", &m_originMomentum, buf_size, split_lvl);
  m_tree->Branch("vertex", &m_vertex, buf_size, split_lvl);
  m_tree->Branch("onSurfaceMomentum", &m_onSurfaceMomentum, buf_size,
                 split_lvl);

  // Track ID
  m_tree->Branch("trackId", &m_trackId, buf_size, split_lvl);
  m_tree->Branch("parentTrackId", &m_parentTrackId, buf_size, split_lvl);
  m_tree->Branch("runId", &m_runId, buf_size, split_lvl);

  // Misc
  m_tree->Branch("isSignal", &m_isSignal, "isSignal/I");

  //------------------------------------------------------------------
  // Initialize the data handles
  m_inputClusters.initialize(m_cfg.inputClusters);
}

ProcessCode RootSimClusterWriter::finalize() {
  if (m_file) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootSimClusterWriter::write(const AlgorithmContext& ctx) {
  auto inputClusters = m_inputClusters(ctx);

  ACTS_DEBUG("Received " << inputClusters.size() << " clusters");
  if (inputClusters.empty()) {
    return ProcessCode::SUCCESS;
  }

  std::lock_guard<std::mutex> lock(m_mutex);

  for (const auto& cluster : inputClusters) {
    const auto& clusterSsl = cluster.sourceLink;

    m_geoCenterGlobal = TVector3(clusterSsl.parametersGlob().x(),
                                 clusterSsl.parametersGlob().y(),
                                 clusterSsl.parametersGlob().z());
    m_geoCenterLocal = TVector2(clusterSsl.parametersLoc().x(),
                                clusterSsl.parametersLoc().y());
    m_geoId = clusterSsl.geometryId().sensitive();
    m_eventId = ctx.eventNumber;

    TArrayD data(4);
    for (std::size_t i = 0; i < 4; i++) {
      data[i] = clusterSsl.covariance()(i);
    }
    m_cov.Use(2, 2, data.GetArray());

    std::vector<TVector3> trackHitsGlobal;
    trackHitsGlobal.reserve(cluster.truthHits.size());

    std::vector<TVector2> trackHitsLocal;
    trackHitsLocal.reserve(cluster.truthHits.size());

    std::vector<int> trackId;
    trackId.reserve(cluster.truthHits.size());

    std::vector<int> parentTrackId;
    parentTrackId.reserve(cluster.truthHits.size());

    std::vector<int> runId;
    runId.reserve(cluster.truthHits.size());

    std::vector<int> charge;
    charge.reserve(cluster.truthHits.size());

    std::vector<int> pdgId;
    pdgId.reserve(cluster.truthHits.size());

    std::vector<TLorentzVector> momenta;
    momenta.reserve(cluster.truthHits.size());

    std::vector<TLorentzVector> originMomenta;
    originMomenta.reserve(cluster.truthHits.size());

    std::vector<TVector3> vertices;
    vertices.reserve(cluster.truthHits.size());

    std::vector<TLorentzVector> onSurfaceMomenta;
    onSurfaceMomenta.reserve(cluster.truthHits.size());

    for (const auto& hit : cluster.truthHits) {
      trackHitsLocal.push_back(TVector2(hit.truthParameters[Acts::eBoundLoc0],
                                        hit.truthParameters[Acts::eBoundLoc1]));

      Acts::Vector3 trackHitGlobal = hit.globalPosition;

      trackHitsGlobal.push_back(
          TVector3(trackHitGlobal.x(), trackHitGlobal.y(), trackHitGlobal.z()));

      trackId.push_back(hit.trackId);
      parentTrackId.push_back(hit.parentTrackId);
      runId.push_back(hit.runId);

      charge.push_back(hit.ipParameters.charge());
      pdgId.push_back(hit.ipParameters.particleHypothesis().absolutePdg());

      TLorentzVector originMom;
      originMom.SetPxPyPzE(
          hit.ipParameters.momentum().x(), hit.ipParameters.momentum().y(),
          hit.ipParameters.momentum().z(),
          std::hypot(hit.ipParameters.absoluteMomentum(),
                     hit.ipParameters.particleHypothesis().mass()));
      originMomenta.push_back(originMom);

      TVector3 vertex(hit.ipParameters.position().x(),
                      hit.ipParameters.position().y(),
                      hit.ipParameters.position().z());
      vertices.push_back(vertex);

      TLorentzVector onSurfaceMom;
      double onSurfP = std::abs(1. / hit.truthParameters[Acts::eBoundQOverP]);
      onSurfaceMom.SetPxPyPzE(
          onSurfP * std::sin(hit.truthParameters[Acts::eBoundTheta]) *
              std::cos(hit.truthParameters[Acts::eBoundPhi]),
          onSurfP * std::sin(hit.truthParameters[Acts::eBoundTheta]) *
              std::sin(hit.truthParameters[Acts::eBoundPhi]),
          onSurfP * std::cos(hit.truthParameters[Acts::eBoundTheta]),
          std::hypot(onSurfP, hit.ipParameters.particleHypothesis().mass()));
      onSurfaceMomenta.push_back(onSurfaceMom);
    }

    m_trackHitsGlobal = trackHitsGlobal;
    m_trackHitsLocal = trackHitsLocal;

    m_trackId = trackId;
    m_parentTrackId = parentTrackId;
    m_runId = runId;

    m_charge = charge;
    m_pdgId = pdgId;

    m_originMomentum = originMomenta;
    m_vertex = vertices;
    m_onSurfaceMomentum = onSurfaceMomenta;

    m_isSignal = cluster.isSignal;

    // Fill the tree
    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
