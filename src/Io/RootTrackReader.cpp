#include "TrackingPipeline/Io/RootTrackReader.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include <Acts/Definitions/PdgParticle.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <stdexcept>
#include <vector>

#include <TFile.h>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

using namespace Acts::UnitLiterals;

RootTrackReader::RootTrackReader(const Config& config,
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

  // m_chain = new TChain(m_cfg.treeName.c_str());
  m_file = new TFile(m_cfg.filePaths.at(0).c_str());
  m_chain = m_file->Get<TTree>(m_cfg.treeName.c_str());

  //------------------------------------------------------------------
  // Tree branches

  // Measurement hits
  m_chain->SetBranchAddress("trackHitsGlobal", &m_trackHitsGlobal);
  m_chain->SetBranchAddress("trackHitsLocal", &m_trackHitsLocal);

  // Covariances of the track hits
  m_chain->SetBranchAddress("trackHitsCovs", &m_trackHitCovs);

  // Geometry ids of the track hits
  m_chain->SetBranchAddress("geometryIds", &m_geometryIds);

  // KF predicted track hits
  m_chain->SetBranchAddress("predictedTrackHitsGlobal",
                            &m_predictedTrackHitsGlobal);
  m_chain->SetBranchAddress("filteredTrackHitsGlobal",
                            &m_filteredTrackHitsGlobal);
  m_chain->SetBranchAddress("smoothedTrackHitsGlobal",
                            &m_smoothedTrackHitsGlobal);

  m_chain->SetBranchAddress("predictedTrackHitsLocal",
                            &m_predictedTrackHitsLocal);
  m_chain->SetBranchAddress("filteredTrackHitsLocal",
                            &m_filteredTrackHitsLocal);
  m_chain->SetBranchAddress("smoothedTrackHitsLocal",
                            &m_smoothedTrackHitsLocal);

  // KF residuals with respect to the measurements
  m_chain->SetBranchAddress("predictedResiduals", &m_predictedResiduals);
  m_chain->SetBranchAddress("filteredResiduals", &m_filteredResiduals);
  m_chain->SetBranchAddress("smoothedResiduals", &m_smoothedResiduals);

  // KF pulls with respect to the measurements
  m_chain->SetBranchAddress("predictedPulls", &m_predictedPulls);
  m_chain->SetBranchAddress("filteredPulls", &m_filteredPulls);
  m_chain->SetBranchAddress("smoothedPulls", &m_smoothedPulls);

  /// Guessed bound track parameters
  m_chain->SetBranchAddress("boundTrackParametersGuess",
                            &m_boundTrackParametersGuess);
  m_chain->SetBranchAddress("boundTrackCovGuess", &m_boundTrackCovGuess);

  /// KF predicted bound track parameters
  m_chain->SetBranchAddress("boundTrackParametersEst",
                            &m_boundTrackParametersEst);
  m_chain->SetBranchAddress("boundTrackCovEst", &m_boundTrackCovEst);

  /// Initial guess of the momentum at the IP
  m_chain->SetBranchAddress("ipMomentumGuess", &m_ipMomentumGuess);

  /// Initial guess of the vertex at the IP
  m_chain->SetBranchAddress("vertexGuess", &m_vertexGuess);

  /// KF predicted momentum at the IP
  m_chain->SetBranchAddress("ipMomentumEst", &m_ipMomentumEst);

  /// KF predicted vertex at the IP
  m_chain->SetBranchAddress("vertexEst", &m_vertexEst);

  // Chi2 and ndf of the fitted track
  m_chain->SetBranchAddress("chi2Predicted", &m_chi2Predicted);
  m_chain->SetBranchAddress("chi2Filtered", &m_chi2Filtered);
  m_chain->SetBranchAddress("chi2Smoothed", &m_chi2Smoothed);
  m_chain->SetBranchAddress("ndf", &m_ndf);

  // Track ID
  m_chain->SetBranchAddress("trackId", &m_trackId);

  // Event ID
  m_chain->SetBranchAddress("eventId", &m_eventId);

  // PDG ID
  m_chain->SetBranchAddress("pdgId", &m_pdgId);

  // Charge
  m_chain->SetBranchAddress("charge", &m_charge);

  // Add the files to the chain
  // for (const auto& path : m_cfg.filePaths) {
  //   m_chain->Add(path.c_str());
  // }

  // Disable all branches and only enable event-id for a first scan of the
  // file
  m_chain->SetBranchStatus("*", false);
  if (!m_chain->GetBranch("eventId")) {
    throw std::invalid_argument("Missing eventId SetbranchAddress");
  }
  m_chain->SetBranchStatus("eventId", true);
  auto nEntries = static_cast<std::size_t>(m_chain->GetEntries());

  // Go through all entries and store the position of the events
  m_chain->GetEntry(0);
  m_eventMap.emplace_back(m_eventId, 0, 0);
  if (!m_cfg.mergeIntoOneEvent) {
    for (std::size_t i = 0; i < nEntries; ++i) {
      m_chain->GetEntry(i);
      if (m_eventId != std::get<0>(m_eventMap.back())) {
        std::get<2>(m_eventMap.back()) = i;
        m_eventMap.emplace_back(m_eventId, i, i);
      }
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
  m_outputSourceLinks.initialize(m_cfg.outputMeasurements);
  m_outputSeedsGuess.initialize(m_cfg.outputSeedsGuess);
  m_outputSeedsEst.initialize(m_cfg.outputSeedsEst);
}

std::pair<std::size_t, std::size_t> RootTrackReader::availableEvents() const {
  return {std::get<0>(m_eventMap.front()), std::get<0>(m_eventMap.back()) + 1};
}

ProcessCode RootTrackReader::read(const AlgorithmContext& ctx) {
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
  Seeds seedsGuess{};
  Seeds seedsEst{};
  std::size_t eventId = std::get<0>(*it);
  std::size_t sslIdx = 0;
  for (auto entry = std::get<1>(*it); entry < std::get<2>(*it); entry++) {
    m_chain->GetEntry(entry);

    // if (m_chi2Smoothed < m_cfg.minChi2 || m_chi2Smoothed > m_cfg.maxChi2) {
    //   continue;
    // }
    // if (m_vertexEst->Y() < -20_mm || m_vertexEst->Y() > 30_mm) {
    //   continue;
    // }
    // if (m_vertexEst->Z() < -40_mm || m_vertexEst->Z() > 40_mm) {
    //   continue;
    // }
    // if (m_ipMomentumEst->P() < 1.5_GeV || m_ipMomentumEst->P() > 3_GeV) {
    //   continue;
    // }

    Acts::BoundMatrix ipCovGuess;
    Acts::BoundMatrix ipCovEst;
    for (std::size_t i = 0; i < Acts::eBoundSize; i++) {
      for (std::size_t j = 0; j < Acts::eBoundSize; j++) {
        ipCovGuess(i, j) = (*m_boundTrackCovGuess)(i, j);
        ipCovEst(i, j) = (*m_boundTrackCovEst)(i, j);
      }
    }

    // Initial guess
    Acts::Vector4 vertexGuess(m_vertexGuess->X(), m_vertexGuess->Y(),
                              m_vertexGuess->Z(), 0);
    Acts::Vector3 ipDirectionGuess(
        m_ipMomentumGuess->X(), m_ipMomentumGuess->Y(), m_ipMomentumGuess->Z());
    ipDirectionGuess.normalize();
    Acts::ParticleHypothesis hypothesis =
        Acts::ParticleHypothesis(Acts::PdgParticle(m_pdgId));
    Acts::CurvilinearTrackParameters ipParametersGuess(
        vertexGuess, ipDirectionGuess, m_charge / m_ipMomentumGuess->P(),
        ipCovGuess, hypothesis);

    // Estimated
    Acts::Vector4 vertexEst(m_vertexEst->X(), m_vertexEst->Y(),
                            m_vertexEst->Z(), 0);
    Acts::Vector3 ipDirectionEst(m_ipMomentumEst->X(), m_ipMomentumEst->Y(),
                                 m_ipMomentumEst->Z());
    ipDirectionEst.normalize();
    Acts::CurvilinearTrackParameters ipParametersEst(
        vertexEst, ipDirectionEst, m_charge / m_ipMomentumEst->P(), ipCovEst,
        hypothesis);

    std::vector<Acts::SourceLink> trackSourceLinks{};
    for (std::size_t i = 0; i < m_trackHitsGlobal->size(); i++) {
      Acts::Vector2 trackHitLocal(m_trackHitsLocal->at(i).X(),
                                  m_trackHitsLocal->at(i).Y());
      Acts::Vector3 trackHitGlobal(m_trackHitsGlobal->at(i).X(),
                                   m_trackHitsGlobal->at(i).Y(),
                                   m_trackHitsGlobal->at(i).Z());
      Acts::SquareMatrix2 trackHitCov;
      trackHitCov << m_trackHitCovs->at(i)(0, 0), m_trackHitCovs->at(i)(0, 1),
          m_trackHitCovs->at(i)(1, 0), m_trackHitCovs->at(i)(1, 1);

      Acts::GeometryIdentifier geoId;
      geoId.setSensitive(m_geometryIds->at(i));
      SimpleSourceLink obsSourceLink(trackHitLocal, trackHitGlobal, trackHitCov,
                                     geoId, eventId, sslIdx);
      sourceLinks.push_back(Acts::SourceLink{obsSourceLink});
      trackSourceLinks.push_back(Acts::SourceLink{obsSourceLink});

      sslIdx++;
    }
    seedsGuess.emplace_back(trackSourceLinks, ipParametersGuess,
                            static_cast<int>(seedsGuess.size()));
    seedsEst.emplace_back(trackSourceLinks, ipParametersEst,
                          static_cast<int>(seedsEst.size()));
  }

  ACTS_DEBUG("Read " << sourceLinks.size() << " source links");
  ACTS_DEBUG("Read " << seedsGuess.size() << " seeds");
  m_outputSourceLinks(ctx, std::move(sourceLinks));
  m_outputSeedsGuess(ctx, std::move(seedsGuess));
  m_outputSeedsEst(ctx, std::move(seedsEst));

  // Return success flag
  return ProcessCode::SUCCESS;
}
