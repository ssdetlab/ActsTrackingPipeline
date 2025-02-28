#include "TrackingPipeline/Io/E320RootSimDataReader.hpp"

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

// Global to local conversion
//
// Here the conventions of all the tranformations
// involved are preserved
Acts::Vector2 convertToLoc(const Acts::Vector3& glob,
                           const Acts::GeometryIdentifier geoId,
                           const E320Geometry::GeometryOptions& gOpt) {
  int nChip = geoId.sensitive() % 10 - 1;
  Acts::Vector2 loc =
      Acts::Vector2(glob.y() - gOpt.chipY.at(nChip), -(glob.x() - gOpt.chipX));
  return loc;
}

inline void E320Io::E320RootSimDataReader::prepareMeasurements(
    const AlgorithmContext& context, std::vector<Acts::SourceLink>* sourceLinks,
    SimClusters* clusters) const {
  // Check if the event number is correct
  auto eventId = m_intColumns.at("eventId");
  if (eventId != context.eventNumber) {
    return;
  }

  // Columns with the measurable quantities
  std::int32_t geoIdval;
  std::int32_t sizeX;
  std::int32_t sizeY;
  std::int32_t size;
  TVector3* geoCenter;

  // Columns with the truth quantities
  std::vector<std::int32_t>* trackId;
  std::vector<std::int32_t>* parentTrackId;
  std::vector<std::int32_t>* runId;
  std::vector<TVector3>* trueHits;
  std::vector<TVector3>* vertices;
  std::vector<TLorentzVector>* mom;
  std::vector<TLorentzVector>* momIP;
  try {
    //--------------------------------
    // Measurable quantities

    // Geometry ID of the chip where
    // the cluster occured
    geoIdval = m_intColumns.at("geoId");

    // The cluster geometrical center
    // in the global coordinates
    geoCenter = m_vector3Columns.at("rglobal_geo");

    // The cluster size in X
    sizeX = m_intColumns.at("xsize");

    // The cluster size in Y
    sizeY = m_intColumns.at("ysize");

    // Cluster size
    size = m_intColumns.at("size");

    //--------------------------------
    // Truth quantities

    // Track IDs of the particles
    // that created the cluster
    trackId = m_vIntColumns.at("tru_trackId");

    // Parent track IDs of the particles
    // that created the cluster
    parentTrackId = m_vIntColumns.at("tru_parenttrackId");

    // Ptarmigan run ID of the particles
    // that created the cluster
    runId = m_vIntColumns.at("tru_runId");

    // The true hit positions of the particles
    // that created the cluster
    trueHits = m_vVector3Columns.at("tru_hit");

    // The true vertex positions of the particles
    // that created the cluster
    vertices = m_vVector3Columns.at("tru_vertex");

    // The true momenta of the particles
    // that created the cluster
    mom = m_vLorentzColumns.at("tru_p");

    // The true momenta of the particles
    // that created the cluster at the IP
    momIP = m_vLorentzColumns.at("tru_p_ip");
  } catch (const std::out_of_range& e) {
    throw std::runtime_error("Missing columns in the ROOT file");
  }

  //-------------------------------
  // Measurable quantities

  // Apply the Geometry ID convention
  Acts::GeometryIdentifier geoId;
  geoId.setSensitive(geoIdval + 11);

  // Create IP covariance matrix from
  // reasonable standard deviations
  Acts::BoundVector ipStdDev;
  ipStdDev[Acts::eBoundLoc0] = 100_um;
  ipStdDev[Acts::eBoundLoc1] = 100_um;
  ipStdDev[Acts::eBoundTime] = 25_ns;
  ipStdDev[Acts::eBoundPhi] = 2_degree;
  ipStdDev[Acts::eBoundTheta] = 2_degree;
  ipStdDev[Acts::eBoundQOverP] = 1 / 100_GeV;
  Acts::BoundSquareMatrix ipCov = ipStdDev.cwiseProduct(ipStdDev).asDiagonal();

  // Convert the cluster geometrical center
  // to the local coordinates
  Acts::Vector3 hitGlob = {geoCenter->X() * Acts::UnitConstants::mm,
                           geoCenter->Y() * Acts::UnitConstants::mm,
                           geoCenter->Z() * Acts::UnitConstants::mm};

  const Acts::Vector2 hitLoc = convertToLoc(hitGlob, geoId, m_cfg.gOpt);

  // Estimate error from the cluster size
  double errX = m_pixSizeX / std::sqrt(12 * size);
  double errY = m_pixSizeY / std::sqrt(12 * size);
  Acts::Vector2 stdDev(errX, errY);
  Acts::SquareMatrix2 cov = stdDev.cwiseProduct(stdDev).asDiagonal();

  // Fill the measurement
  SimpleSourceLink ssl(hitLoc, cov, geoId, eventId, sourceLinks->size());

  //-------------------------------
  // Truth quantities

  // Fill the cluster data
  SimCluster cluster{.sourceLink = ssl};

  if (m_cfg.clusterFilter != nullptr &&
      !m_cfg.clusterFilter->operator()(context.geoContext, cluster)) {
    return;
  }

  sourceLinks->push_back(Acts::SourceLink(ssl));
  double me = 0.511 * Acts::UnitConstants::MeV;
  for (int idx = 0; idx < trueHits->size(); idx++) {
    auto hitMom = mom->at(idx);

    // Convert the vertex
    Acts::Vector3 trueVertex3 = {
        vertices->at(idx).X() * Acts::UnitConstants::mm,
        vertices->at(idx).Y() * Acts::UnitConstants::mm,
        vertices->at(idx).Z() * Acts::UnitConstants::mm};
    trueVertex3 = m_actsToWorld * trueVertex3;

    // KF accepts 4D vectors
    Acts::Vector4 trueVertex = {trueVertex3.x(), trueVertex3.y(),
                                trueVertex3.z(), 0};

    // Convert the true hit to the local coordinates
    Acts::Vector3 trueHitGlob = {
        trueHits->at(idx).X() * Acts::UnitConstants::mm,
        trueHits->at(idx).Y() * Acts::UnitConstants::mm,
        trueHits->at(idx).Z() * Acts::UnitConstants::mm};

    const Acts::Vector2 trueHitLoc =
        convertToLoc(trueHitGlob, geoId, m_cfg.gOpt);

    // Set up the truth parameters
    Acts::BoundVector truthPars = Acts::BoundVector::Zero();
    truthPars[Acts::eBoundLoc0] = trueHitLoc[Acts::eBoundLoc0];
    truthPars[Acts::eBoundLoc1] = trueHitLoc[Acts::eBoundLoc1];

    // Get the momentum at the first hit
    Acts::Vector3 trueP = {hitMom.Px(), hitMom.Py(), hitMom.Pz()};

    Acts::Vector3 dir = trueP.normalized();

    Acts::Vector3 dirRotated = m_actsToWorld * dir;

    truthPars[Acts::eBoundPhi] = Acts::VectorHelpers::phi(dirRotated);
    truthPars[Acts::eBoundTheta] = Acts::VectorHelpers::theta(dirRotated);
    truthPars[Acts::eBoundQOverP] =
        1_e / (hitMom.P() * Acts::UnitConstants::GeV);
    truthPars[Acts::eBoundTime] = hitMom.T();

    // Set up IP information
    TLorentzVector ipMom = momIP->at(idx);

    Acts::Vector3 ipMom3 = {ipMom.Px(), ipMom.Py(), ipMom.Pz()};
    Acts::Vector3 dirIP = ipMom3.normalized();

    Acts::Vector3 dirIPRotated = m_actsToWorld * dirIP;

    Acts::CurvilinearTrackParameters ipParameters(
        trueVertex, Acts::VectorHelpers::phi(dirIPRotated),
        Acts::VectorHelpers::theta(dirIPRotated),
        1_e / (ipMom.P() * Acts::UnitConstants::GeV), ipCov,
        Acts::ParticleHypothesis::electron());

    SimpleSourceLink truthSsl(trueHitLoc, cov, geoId, eventId,
                              clusters->size());

    // Fill the cluster data
    cluster.truthHits.push_back(SimHit(Acts::SourceLink(truthSsl), truthPars,
                                       ipParameters, trackId->at(idx),
                                       parentTrackId->at(idx), runId->at(idx)));
    cluster.isSignal = m_intColumns.at("isSignal");
  }

  clusters->push_back(cluster);
};
