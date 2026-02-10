#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Infrastructure/WriterRegistry.hpp"

#include <toml.hpp>

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
  int bufSize = 32000;
  int splitLvl = 0;

  // Cluster parameters
  m_tree->Branch("geoCenterGlobal", &m_geoCenterGlobal, bufSize, splitLvl);
  m_tree->Branch("geoCenterLocal", &m_geoCenterLocal, bufSize, splitLvl);
  m_tree->Branch("clusterCov", &m_clusterCov, bufSize, splitLvl);
  m_tree->Branch("geoId", &m_geoId, bufSize, splitLvl);

  // Measurement hits
  m_tree->Branch("trackHitsGlobal", &m_trackHitsGlobal, bufSize, splitLvl);
  m_tree->Branch("trackHitsLocal", &m_trackHitsLocal, bufSize, splitLvl);
  m_tree->Branch("eventId", &m_eventId, bufSize, splitLvl);
  m_tree->Branch("charge", &m_charge, bufSize, splitLvl);
  m_tree->Branch("pdgId", &m_pdgId, bufSize, splitLvl);

  // Bound origin parameters
  m_tree->Branch("boundTrackParameters", &m_boundTrackParameters, bufSize,
                 splitLvl);
  m_tree->Branch("boundTrackCov", &m_boundTrackCov, bufSize, splitLvl);

  // Origin momentum
  m_tree->Branch("originMomentum", &m_originMomentum, bufSize, splitLvl);

  // Origin vertex
  m_tree->Branch("vertex", &m_vertex, bufSize, splitLvl);

  // Momentum at clusters
  m_tree->Branch("onSurfaceMomentum", &m_onSurfaceMomentum, bufSize, splitLvl);

  // Track ID
  m_tree->Branch("trackId", &m_trackId, bufSize, splitLvl);
  m_tree->Branch("parentTrackId", &m_parentTrackId, bufSize, splitLvl);
  m_tree->Branch("runId", &m_runId, bufSize, splitLvl);

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

    const Acts::Vector3& clusterParsGlob = clusterSsl.parametersGlob();
    const Acts::Vector2& clusterParsLoc = clusterSsl.parametersLoc();
    m_geoCenterGlobal =
        TVector3(clusterParsGlob.x(), clusterParsGlob.y(), clusterParsGlob.z());
    m_geoCenterLocal = TVector2(clusterParsLoc.x(), clusterParsLoc.y());
    m_geoId = clusterSsl.geometryId().sensitive();
    m_eventId = ctx.eventNumber;

    m_isSignal = cluster.isSignal;

    TArrayD clusterCovData(4);
    for (std::size_t i = 0; i < 4; i++) {
      clusterCovData[i] = clusterSsl.covariance()(i);
    }
    m_clusterCov.Use(2, 2, clusterCovData.GetArray());

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

    m_onSurfaceMomentum.clear();
    m_onSurfaceMomentum.reserve(truthHitsSize);

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
      m_onSurfaceMomentum.push_back(onSurfaceMom);
    }

    // Fill the tree
    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}

namespace {

struct RootSimClusterWriterRegistrar {
  RootSimClusterWriterRegistrar() {
    using namespace TrackingPipeline;

    WriterRegistry::instance().registerBuilder(
      "RootSimClusterWriter",
      [](const toml::value& section,
         Acts::Logging::Level logLevel,
         const std::string& runRoot) -> WriterPtr {

        RootSimClusterWriter::Config cfg;

        cfg.inputClusters =
            toml::find<std::string>(section, "inputClusters");
        cfg.treeName =
            toml::find<std::string>(section, "treeName");
        cfg.filePath =
            runRoot + "/" + toml::find<std::string>(section, "filePath");
            
        return std::make_shared<RootSimClusterWriter>(cfg, logLevel);
      });
  }
} _RootSimClusterWriterRegistrar;

}  // namespace