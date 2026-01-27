#include "TrackingPipeline/Io/ApollonRootSimDataReader.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/ApollonGeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Infrastructure/ReaderRegistry.hpp"

#include <toml.hpp>

namespace ApollonIo {

using go = ApollonGeometry::GeometryOptions;

using namespace Acts::UnitLiterals;

ApollonRootSimDataReader::ApollonRootSimDataReader(const Config& config,
                                                   Acts::Logging::Level level)
    : IReader(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  m_chain = new TChain(m_cfg.treeName.c_str());

  if (m_cfg.filePaths.empty()) {
    throw std::invalid_argument("Missing input filenames");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
  m_outputSimClusters.initialize(m_cfg.outputSimClusters);

  m_chain->SetBranchAddress("geoId", &m_geoId);
  m_chain->SetBranchAddress("pixels", &m_pixels);

  m_chain->SetBranchAddress("isSignal", &m_isSignal);

  m_chain->SetBranchAddress("geoCenterLocal", &m_geoCenterLocal);
  m_chain->SetBranchAddress("geoCenterGlobal", &m_geoCenterGlobal);

  m_chain->SetBranchAddress("totEDep", &m_totEDep);

  m_chain->SetBranchAddress("trackId", &m_trackId);
  m_chain->SetBranchAddress("parentTrackId", &m_parentTrackId);
  std::string idxBranch;
  if (m_cfg.eventSplit) {
    m_chain->SetBranchAddress("eventId", &m_eventId);
    m_chain->SetBranchAddress("runId", &m_runId);
    idxBranch = "eventId";
  } else {
    m_chain->SetBranchAddress("eventId", &m_runId);
    m_chain->SetBranchAddress("runId", &m_eventId);
    idxBranch = "runId";
  }

  m_chain->SetBranchAddress("hitPosGlobal", &m_hitPosGlobal);
  m_chain->SetBranchAddress("hitPosLocal", &m_hitPosLocal);

  m_chain->SetBranchAddress("hitMomDir", &m_hitMomDir);
  m_chain->SetBranchAddress("hitE", &m_hitE);
  m_chain->SetBranchAddress("hitP", &m_hitP);

  m_chain->SetBranchAddress("ipMomDir", &m_ipMomDir);
  m_chain->SetBranchAddress("ipE", &m_ipE);
  m_chain->SetBranchAddress("ipP", &m_ipP);
  m_chain->SetBranchAddress("vertex", &m_vertex);

  m_chain->SetBranchAddress("eDep", &m_eDep);
  m_chain->SetBranchAddress("pdgId", &m_pdgId);

  // Add the files to the chain
  for (const auto& path : m_cfg.filePaths) {
    m_chain->Add(path.c_str());
  }

  // Disable all branches and only enable event-id for a first scan of the
  // file
  m_chain->SetBranchStatus("*", false);
  if (!m_chain->GetBranch(idxBranch.c_str())) {
    throw std::invalid_argument("Missing eventId branch");
  }
  m_chain->SetBranchStatus(idxBranch.c_str(), true);

  auto nEntries = static_cast<std::size_t>(m_chain->GetEntries());

  // Add the first entry
  m_chain->GetEntry(0);
  m_eventMap.emplace_back(m_eventId, 0, 0);

  // Go through all entries and store the position of the events
  for (std::size_t i = 1; i < nEntries; ++i) {
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

  const auto& goInst = *go::instance();
  m_setupRotation =
      Acts::AngleAxis3(goInst.setupRotationAngle,
                       detail::binningValueToDirection(goInst.longBinValue))
          .toRotationMatrix();

  m_vertexCorrection[goInst.primaryIdx] = 0;
  m_vertexCorrection[goInst.longIdx] = goInst.setupLongTranslation;
  m_vertexCorrection[goInst.shortIdx] = 0;

  Acts::BoundVector ipStdDev;
  ipStdDev[Acts::eBoundLoc0] = 1_mm;
  ipStdDev[Acts::eBoundLoc1] = 1_mm;
  ipStdDev[Acts::eBoundTime] = 25_ns;
  ipStdDev[Acts::eBoundPhi] = 2_degree;
  ipStdDev[Acts::eBoundTheta] = 2_degree;
  ipStdDev[Acts::eBoundQOverP] = 1 / 100_GeV;
  m_ipCov = ipStdDev.cwiseProduct(ipStdDev).asDiagonal();

  double errX = 5_um;
  double errY = 5_um;
  Acts::Vector2 stdDev(errX, errY);
  m_hitCov = stdDev.cwiseProduct(stdDev).asDiagonal();
}

std::pair<std::size_t, std::size_t> ApollonRootSimDataReader::availableEvents()
    const {
  return {std::get<0>(m_eventMap.front()), std::get<0>(m_eventMap.back()) + 1};
}

ProcessCode ApollonRootSimDataReader::read(const AlgorithmContext& context) {
  auto it = std::ranges::find_if(m_eventMap, [&](const auto& a) {
    return std::get<0>(a) == context.eventNumber;
  });

  if (it == m_eventMap.end()) {
    // explicitly warn if it happens for the first or last event as that might
    // indicate a human error
    if ((context.eventNumber == availableEvents().first) &&
        (context.eventNumber == availableEvents().second - 1)) {
      ACTS_WARNING("Reading empty event: " << context.eventNumber);
    } else {
      ACTS_DEBUG("Reading empty event: " << context.eventNumber);
    }

    m_outputSourceLinks(context, {});
    m_outputSimClusters(context, {});

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
  SimClusters clusters{};

  for (auto entry = std::get<1>(*it); entry < std::get<2>(*it); entry++) {
    m_chain->GetEntry(entry);

    if (m_eventId != context.eventNumber) {
      continue;
    }

    Acts::GeometryIdentifier geoId;
    geoId.setSensitive(m_geoId);

    Acts::Vector2 hitLoc(m_geoCenterLocal->X() * Acts::UnitConstants::mm,
                         m_geoCenterLocal->Y() * Acts::UnitConstants::mm);
    // Need high precision so use surfaces for conversion
    Acts::Vector3 hitGlob = m_cfg.surfaceMap.at(geoId)->localToGlobal(
        context.geoContext, hitLoc, Acts::Vector3::UnitX());

    SimpleSourceLink ssl(hitLoc, hitGlob, m_hitCov, geoId, m_eventId,
                         sourceLinks.size());

    sourceLinks.emplace_back(ssl);

    SimCluster cluster{.sourceLink = ssl};
    cluster.isSignal = m_isSignal;
    for (std::size_t i = 0; i < m_hitPosLocal->size(); i++) {
      Acts::Vector3 vertex3(m_vertex->at(i).Z() * Acts::UnitConstants::mm,
                            m_vertex->at(i).Y() * Acts::UnitConstants::mm,
                            m_vertex->at(i).X() * Acts::UnitConstants::mm);
      vertex3 = vertex3 - m_vertexCorrection;
      vertex3 = m_setupRotation * vertex3;
      Acts::Vector4 vertex(vertex3.x(), vertex3.y(), vertex3.z(), 0);

      Acts::Vector3 momDir(m_hitMomDir->at(i).Z(), m_hitMomDir->at(i).Y(),
                           m_hitMomDir->at(i).X());
      momDir = m_setupRotation * momDir;

      Acts::BoundVector truthPars = Acts::BoundVector::Zero();
      truthPars[Acts::eBoundLoc0] = hitLoc[Acts::eBoundLoc0];
      truthPars[Acts::eBoundLoc1] = hitLoc[Acts::eBoundLoc1];
      truthPars[Acts::eBoundPhi] = Acts::VectorHelpers::phi(momDir);
      truthPars[Acts::eBoundTheta] = Acts::VectorHelpers::theta(momDir);
      truthPars[Acts::eBoundQOverP] =
          -1_e / (m_hitP->at(i) * Acts::UnitConstants::MeV);
      truthPars[Acts::eBoundTime] = 0;

      Acts::Vector3 momDirIP(m_ipMomDir->at(i).Z(), m_ipMomDir->at(i).Y(),
                             m_ipMomDir->at(i).X());
      momDirIP = m_setupRotation * momDirIP;

      Acts::CurvilinearTrackParameters ipParameters(
          vertex, Acts::VectorHelpers::phi(momDirIP),
          Acts::VectorHelpers::theta(momDirIP),
          -1_e / (m_ipP->at(i) * Acts::UnitConstants::MeV), m_ipCov,
          Acts::ParticleHypothesis::electron());

      cluster.truthHits.emplace_back(std::move(truthPars), std::move(hitGlob),
                                     std::move(ipParameters), m_trackId->at(i),
                                     m_parentTrackId->at(i), m_runId);
    }
    clusters.push_back(cluster);
  }

  ACTS_DEBUG("Sending " << sourceLinks.size() << " measurements");
  ACTS_DEBUG("Sending " << clusters.size() << " sim clusters");

  m_outputSourceLinks(context, std::move(sourceLinks));
  m_outputSimClusters(context, std::move(clusters));

  // Return success flag
  return ProcessCode::SUCCESS;
}

}  // namespace ApollonIo

namespace {
  struct ApollonRootSimDataReaderRegistrar {
    ApollonRootSimDataReaderRegistrar() {
      using namespace TrackingPipeline;

      ReaderRegistry::instance().registerBuilder(
          "ApollonRootSimDataReader",
          [](const toml::value& section,
             const SurfaceMap& surfaceMap,
             Acts::Logging::Level logLevel) -> ReaderPtr {

            ApollonIo::ApollonRootSimDataReader::Config cfg;

            // Fill config from TOML section
            cfg.outputSourceLinks =
                toml::find<std::string>(section, "outputSourceLinks");
            cfg.outputSimClusters =
                toml::find<std::string>(section, "outputSimClusters");
            cfg.filePaths =
                toml::find<std::vector<std::string>>(section, "filePaths");
            cfg.treeName =
                toml::find<std::string>(section, "treeName");
            cfg.eventSplit =
                toml::find<bool>(section, "eventSplit", false);
            cfg.surfaceMap = surfaceMap;

            return std::make_shared<ApollonIo::ApollonRootSimDataReader>(cfg, logLevel);
          });
    }
  } _ApollonRootSimDataReaderRegistrar;
} //namespace