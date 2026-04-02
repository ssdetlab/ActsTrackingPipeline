#include "TrackingPipeline/Io/RootSimTrackReader.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include <Acts/Definitions/PdgParticle.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "TFile.h"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

using namespace Acts::UnitLiterals;

RootSimTrackReader::RootSimTrackReader(const Config& config,
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

  // True hits
  m_tree->SetBranchAddress("trueTrackHitsGlobal", &m_trueTrackHitsGlobal);
  m_tree->SetBranchAddress("trueTrackHitsLocal", &m_trueTrackHitsLocal);
  m_tree->SetBranchAddress("onSurfaceMomentumTruth", &m_onSurfaceMomentumTruth);
  m_tree->SetBranchAddress("isSignal", &m_isSignal);

  // Measurement hits
  m_tree->SetBranchAddress("trackHitsGlobal", &m_trackHitsGlobal);
  m_tree->SetBranchAddress("trackHitsLocal", &m_trackHitsLocal);

  // Covariances of the track hits
  m_tree->SetBranchAddress("trackHitsCovs", &m_trackHitCovs);

  // Geometry ids of the track hits
  m_tree->SetBranchAddress("geometryIds", &m_geometryIds);

  // KF predicted track hits
  m_tree->SetBranchAddress("predictedTrackHitsGlobal",
                           &m_predictedTrackHitsGlobal);
  m_tree->SetBranchAddress("filteredTrackHitsGlobal",
                           &m_filteredTrackHitsGlobal);
  m_tree->SetBranchAddress("smoothedTrackHitsGlobal",
                           &m_smoothedTrackHitsGlobal);

  m_tree->SetBranchAddress("predictedTrackHitsLocal",
                           &m_predictedTrackHitsLocal);
  m_tree->SetBranchAddress("filteredTrackHitsLocal", &m_filteredTrackHitsLocal);
  m_tree->SetBranchAddress("smoothedTrackHitsLocal", &m_smoothedTrackHitsLocal);

  // KF residuals with respect to the true hits
  m_tree->SetBranchAddress("truePredictedResiduals", &m_truePredictedResiduals);
  m_tree->SetBranchAddress("trueFilteredResiduals", &m_trueFilteredResiduals);
  m_tree->SetBranchAddress("trueSmoothedResiduals", &m_trueSmoothedResiduals);

  // KF residuals with respect to the measurements
  m_tree->SetBranchAddress("predictedResiduals", &m_predictedResiduals);
  m_tree->SetBranchAddress("filteredResiduals", &m_filteredResiduals);
  m_tree->SetBranchAddress("smoothedResiduals", &m_smoothedResiduals);

  // KF pulls with respect to the true hits
  m_tree->SetBranchAddress("truePredictedPulls", &m_truePredictedPulls);
  m_tree->SetBranchAddress("trueFilteredPulls", &m_trueFilteredPulls);
  m_tree->SetBranchAddress("trueSmoothedPulls", &m_trueSmoothedPulls);

  // KF pulls with respect to the measurements
  m_tree->SetBranchAddress("predictedPulls", &m_predictedPulls);
  m_tree->SetBranchAddress("filteredPulls", &m_filteredPulls);
  m_tree->SetBranchAddress("smoothedPulls", &m_smoothedPulls);

  /// Guessed bound track parameters
  m_tree->SetBranchAddress("boundTrackParametersGuess",
                           &m_boundTrackParametersGuess);
  m_tree->SetBranchAddress("boundTrackCovGuess", &m_boundTrackCovGuess);

  /// KF predicted bound track parameters
  m_tree->SetBranchAddress("boundTrackParametersEst",
                           &m_boundTrackParametersEst);
  m_tree->SetBranchAddress("boundTrackCovEst", &m_boundTrackCovEst);

  /// True bound track parameters
  m_tree->SetBranchAddress("boundTrackParametersTruth",
                           &m_boundTrackParametersTruth);
  m_tree->SetBranchAddress("boundTrackCovTruth", &m_boundTrackCovTruth);

  /// Initial guess of the momentum at the IP
  m_tree->SetBranchAddress("ipMomentumGuess", &m_ipMomentumGuess);

  /// Initial guess of the vertex at the IP
  m_tree->SetBranchAddress("vertexGuess", &m_vertexGuess);

  /// KF predicted momentum at the IP
  m_tree->SetBranchAddress("ipMomentumEst", &m_ipMomentumEst);

  /// KF predicted vertex at the IP
  m_tree->SetBranchAddress("vertexEst", &m_vertexEst);

  // True momentum at the IP
  m_tree->SetBranchAddress("ipMomentumTruth", &m_ipMomentumTruth);

  // True vertex at the IP
  m_tree->SetBranchAddress("vertexTruth", &m_vertexTruth);

  // Chi2 and ndf of the fitted track
  m_tree->SetBranchAddress("chi2Predicted", &m_chi2Predicted);
  m_tree->SetBranchAddress("chi2Filtered", &m_chi2Filtered);
  m_tree->SetBranchAddress("chi2Smoothed", &m_chi2Smoothed);
  m_tree->SetBranchAddress("ndf", &m_ndf);

  // Track ID
  m_tree->SetBranchAddress("stateTrackId", &m_stateTrackId);
  m_tree->SetBranchAddress("stateParentTrackId", &m_stateParentTrackId);
  m_tree->SetBranchAddress("stateRunId", &m_stateRunId);

  m_tree->SetBranchAddress("trackId", &m_trackId);
  m_tree->SetBranchAddress("parentTrackId", &m_parentTrackId);
  m_tree->SetBranchAddress("runId", &m_runId);

  // Event ID
  m_tree->SetBranchAddress("eventId", &m_eventId);

  // True track size
  m_tree->SetBranchAddress("trueTrackSize", &m_trueTrackSize);
  m_tree->SetBranchAddress("capturedTrackSize", &m_capturedTrackSize);

  // PDG ID
  m_tree->SetBranchAddress("pdgId", &m_pdgId);

  // Charge
  m_tree->SetBranchAddress("charge", &m_charge);

  // Disable all branches and only enable event-id for a first scan of the
  // file
  m_tree->SetBranchStatus("*", false);
  if (m_tree->GetBranch("eventId") == nullptr) {
    throw std::invalid_argument("Missing eventId SetbranchAddress");
  }
  m_tree->SetBranchStatus("eventId", true);
  auto nEntries = static_cast<std::size_t>(m_tree->GetEntries());

  // Go through all entries and store the position of the events
  m_tree->GetEntry(0);
  m_eventMap.emplace_back(m_eventId, 0, 0);
  if (!m_cfg.mergeIntoOneEvent) {
    for (std::size_t i = 0; i < nEntries; ++i) {
      m_tree->GetEntry(i);
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
  m_tree->SetBranchStatus("*", true);
  ACTS_DEBUG("Event range: " << availableEvents().first << " - "
                             << availableEvents().second);

  // Initialize the data handles
  m_outputSourceLinks.initialize(m_cfg.outputMeasurements);
  m_outputSimClusters.initialize(m_cfg.outputSimClusters);
  m_outputSeedsGuess.initialize(m_cfg.outputSeedsGuess);
  m_outputSeedsEst.initialize(m_cfg.outputSeedsEst);
}

std::pair<std::size_t, std::size_t> RootSimTrackReader::availableEvents()
    const {
  return {std::get<0>(m_eventMap.front()), std::get<0>(m_eventMap.back()) + 1};
}

ProcessCode RootSimTrackReader::read(const AlgorithmContext& ctx) {
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
  Seeds seedsGuess{};
  Seeds seedsEst{};
  std::size_t eventId = std::get<0>(*it);
  std::size_t sslIdx = 0;
  const Constraints& constraints = m_cfg.constraints;
  for (auto entry = std::get<1>(*it); entry < std::get<2>(*it); entry++) {
    m_tree->GetEntry(entry);

    if (m_chi2Smoothed < constraints.minChi2 ||
        m_chi2Smoothed > constraints.maxChi2) {
      continue;
    }
    if (m_vertexEst->Y() < constraints.minVertexEstLong ||
        m_vertexEst->Y() > constraints.maxVertexEstLong) {
      continue;
    }
    if (m_vertexEst->Z() < constraints.minVertexEstShort ||
        m_vertexEst->Z() > constraints.maxVertexEstShort) {
      continue;
    }
    if (m_ipMomentumEst->P() < constraints.minAbsMomentumEst ||
        m_ipMomentumEst->P() > constraints.maxAbsMomentumEst) {
      continue;
    }

    Acts::BoundMatrix ipCovGuess;
    Acts::BoundMatrix ipCovEst;
    Acts::BoundMatrix ipCovTruth;
    for (std::size_t i = 0; i < Acts::eBoundSize; i++) {
      for (std::size_t j = 0; j < Acts::eBoundSize; j++) {
        ipCovGuess(i, j) = (*m_boundTrackCovGuess)(i, j);
        ipCovEst(i, j) = (*m_boundTrackCovEst)(i, j);
        ipCovTruth(i, j) = (*m_boundTrackCovTruth)(i, j);
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

    // Truth
    Acts::Vector4 vertexTruth(m_vertexTruth->X(), m_vertexTruth->Y(),
                              m_vertexTruth->Z(), 0);
    Acts::Vector3 ipDirectionTruth(
        m_ipMomentumTruth->X(), m_ipMomentumTruth->Y(), m_ipMomentumTruth->Z());
    ipDirectionTruth.normalize();
    Acts::CurvilinearTrackParameters ipParametersTruth(
        vertexTruth, ipDirectionTruth, m_charge / m_ipMomentumTruth->P(),
        ipCovTruth, hypothesis);

    std::vector<Acts::SourceLink> trackSourceLinks{};
    for (std::size_t i = 0; i < m_trueTrackHitsGlobal->size(); i++) {
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

      Acts::BoundVector truthParameters;
      truthParameters[Acts::eBoundLoc0] = m_trueTrackHitsLocal->at(i).X();
      truthParameters[Acts::eBoundLoc1] = m_trueTrackHitsLocal->at(i).Y();
      truthParameters[Acts::eBoundPhi] = m_onSurfaceMomentumTruth->at(i).Phi();
      truthParameters[Acts::eBoundTheta] =
          m_onSurfaceMomentumTruth->at(i).Theta();
      truthParameters[Acts::eBoundQOverP] =
          m_charge / m_onSurfaceMomentumTruth->at(i).P();

      Acts::Vector3 trueTrackHitGlobal(m_trueTrackHitsGlobal->at(i).X(),
                                       m_trueTrackHitsGlobal->at(i).Y(),
                                       m_trueTrackHitsGlobal->at(i).Z());
      SimHit hit{truthParameters,
                 trueTrackHitGlobal,
                 ipParametersTruth,
                 static_cast<int>(m_stateTrackId->at(i)),
                 static_cast<int>(m_stateParentTrackId->at(i)),
                 static_cast<int>(m_stateRunId->at(i))};
      SimCluster cluster{
          obsSourceLink, {hit}, static_cast<bool>(m_isSignal->at(i))};
      simClusters.push_back(cluster);
    }
    seedsGuess.emplace_back(trackSourceLinks, ipParametersGuess,
                            static_cast<int>(seedsGuess.size()));
    seedsEst.emplace_back(trackSourceLinks, ipParametersEst,
                          static_cast<int>(seedsEst.size()));
  }

  ACTS_DEBUG("Read " << sourceLinks.size() << " source links");
  ACTS_DEBUG("Read " << simClusters.size() << " clusters");
  ACTS_DEBUG("Read " << seedsGuess.size() << " seeds");
  m_outputSourceLinks(ctx, std::move(sourceLinks));
  m_outputSimClusters(ctx, std::move(simClusters));
  m_outputSeedsGuess(ctx, std::move(seedsGuess));
  m_outputSeedsEst(ctx, std::move(seedsEst));

  // Return success flag
  return ProcessCode::SUCCESS;
}
