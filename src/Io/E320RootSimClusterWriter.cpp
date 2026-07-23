#include "TrackingPipeline/Io/E320RootSimClusterWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "TrackingPipeline/EventData/ExtendedSourceLink.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

E320::E320RootSimClusterWriter::E320RootSimClusterWriter(
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
  // Tree branches
  int bufSize = 32000;
  int splitLvl = 0;

  // Cluster hit in the surface frame
  m_tree->Branch("geoCenterLocal", &m_geoCenterLocal, bufSize, splitLvl);

  // Cluster hit in the global frame
  m_tree->Branch("geoCenterGlobal", &m_geoCenterGlobal, bufSize, splitLvl);

  // Cluster hit covariance in the surface frame
  m_tree->Branch("clusterCov", &m_clusterCov, bufSize, splitLvl);

  // Momentum direction measurement in the track frame
  m_tree->Branch("onSurfaceDirection", &m_onSurfaceDirection, bufSize,
                 splitLvl);

  // Track angular covariance in the surface frame
  m_tree->Branch("angleCov", &m_angleCov, bufSize, splitLvl);

  // Surface geometry ID
  m_tree->Branch("geoId", &m_geoId, bufSize, splitLvl);

  // Event ID
  m_tree->Branch("eventId", &m_eventId, bufSize, splitLvl);

  // Signal/background flag
  m_tree->Branch("isSignal", &m_isSignal, "isSignal/I");

  // Track IDs
  m_tree->Branch("trackId", &m_trackId, bufSize, splitLvl);
  m_tree->Branch("parentTrackId", &m_parentTrackId, bufSize, splitLvl);
  m_tree->Branch("runId", &m_runId, bufSize, splitLvl);

  // True track hits within the clusters
  m_tree->Branch("trackHitsGlobal", &m_trackHitsGlobal, bufSize, splitLvl);
  m_tree->Branch("trackHitsLocal", &m_trackHitsLocal, bufSize, splitLvl);

  // Bound origin parameters
  m_tree->Branch("boundTrackParameters", &m_boundTrackParameters, bufSize,
                 splitLvl);
  m_tree->Branch("boundTrackCov", &m_boundTrackCov, bufSize, splitLvl);

  // Origin momentum
  m_tree->Branch("originMomentum", &m_originMomentum, bufSize, splitLvl);

  // Origin vertex
  m_tree->Branch("vertex", &m_vertex, bufSize, splitLvl);

  // Momentum at clusters
  m_tree->Branch("onSurfaceMomentumTruth", &m_onSurfaceMomentumTruth, bufSize,
                 splitLvl);

  // Charges of the tracks' particles
  m_tree->Branch("charge", &m_charge, bufSize, splitLvl);

  // PDG IDs of the tracks' particles
  m_tree->Branch("pdgId", &m_pdgId, bufSize, splitLvl);

  //------------------------------------------------------------------
  // Initialize the data handles
  m_inputSimClusters.initialize(m_cfg.inputClusters);
}

ProcessCode E320::E320RootSimClusterWriter::finalize() {
  if (m_file != nullptr) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode E320::E320RootSimClusterWriter::write(const AlgorithmContext& ctx) {
  const auto& inputClusters = m_inputSimClusters(ctx);

  ACTS_DEBUG("Received " << inputClusters.size() << " clusters");
  if (inputClusters.empty()) {
    return ProcessCode::SUCCESS;
  }

  std::lock_guard<std::mutex> lock(m_mutex);

  for (const auto& cluster : inputClusters) {
    const auto& clusterSl = cluster.sourceLink;

    if (clusterSl.type() == typeid(SimpleSourceLink)) {
      const auto& clusterSsl = clusterSl.get<SimpleSourceLink>();
      const Acts::Vector3& clusterParsGlob = clusterSsl.parametersGlob();
      const Acts::Vector2& clusterParsLoc = clusterSsl.parametersLoc();

      m_geoCenterGlobal = TVector3(clusterParsGlob.x(), clusterParsGlob.y(),
                                   clusterParsGlob.z());
      m_geoCenterLocal = TVector2(clusterParsLoc.x(), clusterParsLoc.y());
      m_onSurfaceDirection = TVector3(0, 0, 0);
      m_geoId = clusterSsl.geometryId().sensitive();

      TArrayD clusterCovData(4);
      for (std::size_t i = 0; i < 4; i++) {
        clusterCovData[i] = clusterSsl.covariance()(i);
      }
      m_clusterCov.Use(2, 2, clusterCovData.GetArray());
    } else if (clusterSl.type() == typeid(ExtendedSourceLink)) {
      const auto& clusterEsl = clusterSl.get<ExtendedSourceLink>();
      const Acts::ActsVector<ExtendedGlobalSize>& clusterParsGlob =
          clusterEsl.parametersGlob();
      const Acts::ActsVector<ExtendedLocalSize>& clusterParsLoc =
          clusterEsl.parametersLoc();

      m_geoCenterGlobal =
          TVector3(clusterParsGlob(0), clusterParsGlob(1), clusterParsGlob(2));
      m_geoCenterLocal = TVector2(clusterParsLoc(0), clusterParsLoc(1));
      m_onSurfaceDirection =
          TVector3(clusterParsGlob(3), clusterParsGlob(4), clusterParsGlob(5));
      m_geoId = clusterEsl.geometryId().sensitive();

      Acts::SquareMatrix2 hitCov = clusterEsl.covariance().topLeftCorner(2, 2);
      TArrayD hitCovData(4);
      for (std::size_t i = 0; i < 4; i++) {
        hitCovData[i] = hitCov(i);
      }
      m_clusterCov.Use(2, 2, hitCovData.GetArray());

      Acts::SquareMatrix2 angleCov =
          clusterEsl.covariance().bottomRightCorner(2, 2);
      TArrayD angleCovData(4);
      for (std::size_t i = 0; i < 4; i++) {
        angleCovData[i] = angleCov(i);
      }
      m_angleCov.Use(2, 2, angleCovData.GetArray());
    }

    m_eventId = ctx.eventNumber;

    m_isSignal = static_cast<int>(cluster.isSignal);

    std::size_t truthHitsSize = cluster.truthHits.size();
    m_trackHitsGlobal.clear();
    m_trackHitsGlobal.reserve(truthHitsSize);

    m_trackHitsLocal.clear();
    m_trackHitsLocal.reserve(truthHitsSize);

    m_trackId.clear();
    m_trackId.reserve(truthHitsSize);

    m_parentTrackId.clear();
    m_parentTrackId.reserve(truthHitsSize);

    m_runId.clear();
    m_runId.reserve(truthHitsSize);

    m_charge.clear();
    m_charge.reserve(truthHitsSize);

    m_pdgId.clear();
    m_pdgId.reserve(truthHitsSize);

    m_boundTrackParameters.clear();
    m_boundTrackParameters.reserve(truthHitsSize);

    m_boundTrackCov.clear();
    m_boundTrackCov.reserve(truthHitsSize);

    m_originMomentum.clear();
    m_originMomentum.reserve(truthHitsSize);

    m_vertex.clear();
    m_vertex.reserve(truthHitsSize);

    m_onSurfaceMomentumTruth.clear();
    m_onSurfaceMomentumTruth.reserve(truthHitsSize);

    for (const auto& hit : cluster.truthHits) {
      const Acts::BoundVector& hitTruthPars = hit.truthParameters;
      m_trackHitsLocal.push_back(TVector2(hitTruthPars[Acts::eBoundLoc0],
                                          hitTruthPars[Acts::eBoundLoc1]));

      const Acts::Vector3& trackHitGlobal = hit.globalPosition;

      m_trackHitsGlobal.push_back(
          TVector3(trackHitGlobal.x(), trackHitGlobal.y(), trackHitGlobal.z()));

      m_trackId.push_back(hit.trackId);
      m_parentTrackId.push_back(hit.parentTrackId);
      m_runId.push_back(hit.runId);

      m_charge.push_back(hit.ipParameters.charge());
      m_pdgId.push_back(hit.ipParameters.particleHypothesis().absolutePdg());

      const Acts::CurvilinearTrackParameters& originParameters =
          hit.ipParameters;

      // Origin momentum
      const Acts::Vector3& ipMomentum = originParameters.momentum();
      double particleMass = originParameters.particleHypothesis().mass();

      TLorentzVector originMomentum;
      originMomentum.SetPxPyPzE(ipMomentum.x(), ipMomentum.y(), ipMomentum.z(),
                                std::hypot(ipMomentum.norm(), particleMass));
      m_originMomentum.push_back(originMomentum);

      // Vertex
      const auto& ipPosition = originParameters.position(ctx.geoContext);
      TVector3 vertex(ipPosition.x(), ipPosition.y(), ipPosition.z());
      m_vertex.push_back(vertex);

      // Bound track parameters
      Acts::BoundVector boundTrackParameters = originParameters.parameters();

      TVectorD boundTrackParametersVec;
      TArrayD boundTrackParsData(Acts::eBoundSize);
      for (std::size_t i = 0; i < Acts::eBoundSize; i++) {
        boundTrackParsData[i] = boundTrackParameters(i);
      }
      boundTrackParametersVec.Use(1, Acts::eBoundSize,
                                  boundTrackParsData.GetArray());
      m_boundTrackParameters.push_back(boundTrackParametersVec);

      Acts::BoundMatrix boundTrackCov = originParameters.covariance().value();
      TMatrixD boundTrackCovMat;
      TArrayD boundTrackCovData(Acts::eBoundSize * Acts::eBoundSize);
      for (std::size_t i = 0; i < Acts::eBoundSize * Acts::eBoundSize; i++) {
        boundTrackCovData[i] = boundTrackCov(i);
      }
      boundTrackCovMat.Use(Acts::eBoundSize, Acts::eBoundSize,
                           boundTrackCovData.GetArray());
      m_boundTrackCov.push_back(boundTrackCovMat);

      // On surface momentum
      TLorentzVector onSurfaceMom;
      double onSurfP = std::abs(1. / hitTruthPars[Acts::eBoundQOverP]);
      onSurfaceMom.SetPxPyPzE(
          onSurfP * std::sin(hitTruthPars[Acts::eBoundTheta]) *
              std::cos(hitTruthPars[Acts::eBoundPhi]),
          onSurfP * std::sin(hitTruthPars[Acts::eBoundTheta]) *
              std::sin(hitTruthPars[Acts::eBoundPhi]),
          onSurfP * std::cos(hitTruthPars[Acts::eBoundTheta]),
          std::hypot(onSurfP, hit.ipParameters.particleHypothesis().mass()));
      m_onSurfaceMomentumTruth.push_back(onSurfaceMom);
    }

    // Fill the tree
    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
