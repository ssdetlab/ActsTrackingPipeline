#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/EventData/SourceLink.hpp>

#include <stdexcept>
#include <vector>

#include <TLorentzVector.h>
#include <TVector3.h>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

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
  m_tree->Branch("geoCenter", &m_geoCenter, buf_size, split_lvl);
  m_tree->Branch("trackHits", &m_trackHits, buf_size, split_lvl);
  m_tree->Branch("onSurfMomemtum", &m_onSurfMomemtum, buf_size, split_lvl);

  // Parameters at the origin
  m_tree->Branch("originMomentum", &m_originMomentum, buf_size, split_lvl);
  m_tree->Branch("vertex", &m_vertex, buf_size, split_lvl);

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

  std::lock_guard<std::mutex> lock(m_mutex);

  for (const auto& cluster : inputClusters) {
    const auto& clusterSsl = cluster.sourceLink;

    const auto* surf = m_cfg.surfaceAccessor(Acts::SourceLink(clusterSsl));

    Acts::Vector3 geoCenter = surf->localToGlobal(
        ctx.geoContext, clusterSsl.parameters(), Acts::Vector3(0, 1, 0));

    m_geoCenter = TVector3(geoCenter.x(), geoCenter.y(), geoCenter.z());

    std::vector<TVector3> trackHits;
    trackHits.reserve(cluster.truthHits.size());

    std::vector<int> trackId;
    trackId.reserve(cluster.truthHits.size());

    std::vector<int> parentTrackId;
    parentTrackId.reserve(cluster.truthHits.size());

    std::vector<int> runId;
    runId.reserve(cluster.truthHits.size());

    std::vector<TLorentzVector> momenta;
    momenta.reserve(cluster.truthHits.size());

    std::vector<TLorentzVector> originMomenta;
    originMomenta.reserve(cluster.truthHits.size());

    std::vector<TVector3> vertices;
    vertices.reserve(cluster.truthHits.size());

    for (const auto& hit : cluster.truthHits) {
      const auto& hitSsl = hit.sourceLink.get<SimpleSourceLink>();

      Acts::Vector3 trackHit = surf->localToGlobal(
          ctx.geoContext, hitSsl.parameters(), Acts::Vector3(0, 1, 0));

      trackHits.push_back(TVector3(trackHit.x(), trackHit.y(), trackHit.z()));

      trackId.push_back(hit.trackId);
      parentTrackId.push_back(hit.parentTrackId);
      runId.push_back(hit.runId);

      Acts::FreeVector onSurfParameters = Acts::transformBoundToFreeParameters(
          *surf, ctx.geoContext, hit.truthParameters);

      double onSurfE = std::abs(onSurfParameters[Acts::eFreeQOverP]);

      Acts::Vector3 onSurfMomentum =
          onSurfE * Acts::Vector3(onSurfParameters[Acts::eFreeDir0],
                                  onSurfParameters[Acts::eFreeDir1],
                                  onSurfParameters[Acts::eFreeDir2]);

      TLorentzVector mom;
      mom.SetPxPyPzE(onSurfMomentum.x(), onSurfMomentum.y(), onSurfMomentum.z(),
                     onSurfE);
      momenta.push_back(mom);

      TLorentzVector originMom;
      originMom.SetPxPyPzE(
          hit.ipParameters.momentum().x(), hit.ipParameters.momentum().y(),
          hit.ipParameters.momentum().z(), hit.ipParameters.absoluteMomentum());
      originMomenta.push_back(originMom);

      TVector3 vertex(hit.ipParameters.position().x(),
                      hit.ipParameters.position().y(),
                      hit.ipParameters.position().z());
      vertices.push_back(vertex);
    }

    m_trackHits = trackHits;

    m_trackId = trackId;
    m_parentTrackId = parentTrackId;
    m_runId = runId;

    m_onSurfMomemtum = momenta;
    m_originMomentum = originMomenta;
    m_vertex = vertices;

    m_isSignal = cluster.isSignal;

    // Fill the tree
    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
