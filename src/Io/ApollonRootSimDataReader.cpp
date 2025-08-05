#include "TrackingPipeline/Io/ApollonRootSimDataReader.hpp"

#include <Acts/Definitions/Algebra.hpp>

#include <cstddef>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

inline void ApollonIo::ApollonRootSimDataReader::prepareMeasurements(
    const AlgorithmContext& context, std::vector<Acts::SourceLink>* sourceLinks,
    SimClusters* clusters) const {
  auto eventId = m_intColumns.at("eventId");
  if (eventId != context.eventNumber) {
    return;
  }

  std::vector<int>* geoIdVal;
  std::vector<int>* isSignal;

  std::vector<int>* trackId;
  std::vector<int>* parentTrackId;
  std::vector<int>* runId;

  std::vector<TVector3>* hitPosGlobal;
  std::vector<TVector2>* hitPosLocal;

  std::vector<TVector3>* hitMomDir;
  std::vector<double>* hitE;

  std::vector<TVector3>* ipMomDir;
  std::vector<double>* ipE;
  std::vector<TVector3>* vertices;
  try {
    geoIdVal = m_vIntColumns.at("geoId");
    isSignal = m_vIntColumns.at("isSignal");

    hitPosGlobal = m_vVector3Columns.at("hitPosGlobal");
    hitPosLocal = m_vVector2Columns.at("hitPosLocal");

    trackId = m_vIntColumns.at("trackId");
    parentTrackId = m_vIntColumns.at("parentTrackId");
    runId = m_vIntColumns.at("runId");

    ipMomDir = m_vVector3Columns.at("ipMomDir");
    ipE = m_vDoubleColumns.at("ipE");
    vertices = m_vVector3Columns.at("vertex");

    hitMomDir = m_vVector3Columns.at("hitMomDir");
    hitE = m_vDoubleColumns.at("hitE");
  } catch (const std::out_of_range& /*err*/) {
    throw std::runtime_error("Missing columns in the ROOT file");
  }

  Acts::BoundVector ipStdDev;
  ipStdDev[Acts::eBoundLoc0] = 100_um;
  ipStdDev[Acts::eBoundLoc1] = 100_um;
  ipStdDev[Acts::eBoundTime] = 25_ns;
  ipStdDev[Acts::eBoundPhi] = 2_degree;
  ipStdDev[Acts::eBoundTheta] = 2_degree;
  ipStdDev[Acts::eBoundQOverP] = 1 / 100_GeV;
  Acts::BoundSquareMatrix ipCov = ipStdDev.cwiseProduct(ipStdDev).asDiagonal();

  double errX = 5_um;
  double errY = 5_um;
  Acts::Vector2 stdDev(errX, errY);
  Acts::SquareMatrix2 cov = stdDev.cwiseProduct(stdDev).asDiagonal();

  for (std::size_t i = 0; i < geoIdVal->size(); i++) {
    bool isVcExit = (geoIdVal->at(i) == 100);

    if (isVcExit) {
      continue;
    }

    Acts::GeometryIdentifier geoId;
    geoId.setSensitive(geoIdVal->at(i));

    Acts::Vector2 hitLoc(hitPosLocal->at(i).X(), hitPosLocal->at(i).Y());

    SimpleSourceLink ssl(hitLoc, cov, geoId, eventId, sourceLinks->size());

    std::cout << "-----------------------------------------\n";
    std::cout << "GEO ID " << geoIdVal->at(i) << "\n";
    std::cout << "HIT LOC " << hitLoc.transpose() << "\n";
    std::cout << "HIT GLOB G4 [" << hitPosGlobal->at(i).X() << ", "
              << hitPosGlobal->at(i).Y() << ", " << hitPosGlobal->at(i).Z()
              << "]\n";
    std::cout << "HIT GLOB ACTS "
              << m_cfg.surfaceAccessor(Acts::SourceLink(ssl))
                     ->localToGlobal(context.geoContext, hitLoc,
                                     Acts::Vector3::UnitX())
                     .transpose()
              << "\n";

    SimCluster cluster{.sourceLink = ssl};

    if (m_cfg.clusterFilter != nullptr &&
        !m_cfg.clusterFilter->operator()(context.geoContext, cluster)) {
      return;
    }

    sourceLinks->push_back(Acts::SourceLink(ssl));

    Acts::Vector4 vertex(vertices->at(i).X() * Acts::UnitConstants::mm,
                         vertices->at(i).Y() * Acts::UnitConstants::mm,
                         vertices->at(i).Z() * Acts::UnitConstants::mm, 0);

    Acts::Vector2 trueHitLoc = hitLoc;

    Acts::Vector3 momDir(hitMomDir->at(i).X(), hitMomDir->at(i).Y(),
                         hitMomDir->at(i).Z());

    Acts::BoundVector truthPars = Acts::BoundVector::Zero();
    truthPars[Acts::eBoundLoc0] = trueHitLoc[Acts::eBoundLoc0];
    truthPars[Acts::eBoundLoc1] = trueHitLoc[Acts::eBoundLoc1];
    truthPars[Acts::eBoundPhi] = Acts::VectorHelpers::phi(momDir);
    truthPars[Acts::eBoundTheta] = Acts::VectorHelpers::theta(momDir);
    truthPars[Acts::eBoundQOverP] =
        1_e / (hitE->at(i) * Acts::UnitConstants::MeV);
    truthPars[Acts::eBoundTime] = 0;

    Acts::Vector3 momDirIP(ipMomDir->at(i).X(), ipMomDir->at(i).Y(),
                           ipMomDir->at(i).Z());

    Acts::CurvilinearTrackParameters ipParameters(
        vertex, Acts::VectorHelpers::phi(momDirIP),
        Acts::VectorHelpers::theta(momDirIP),
        1_e / (ipE->at(i) * Acts::UnitConstants::MeV), ipCov,
        Acts::ParticleHypothesis::electron());

    SimpleSourceLink truthSsl(trueHitLoc, cov, geoId, eventId,
                              clusters->size());

    cluster.truthHits.push_back(SimHit(truthPars, ipParameters, trackId->at(i),
                                       parentTrackId->at(i), runId->at(i)));
    cluster.isSignal = isSignal->at(i);
    clusters->push_back(cluster);
  }
}
