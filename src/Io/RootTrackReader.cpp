#include "TrackingPipeline/Io/RootTrackReader.hpp"

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

  // KF residuals with respect to the measurements
  m_tree->SetBranchAddress("predictedResiduals", &m_predictedResiduals);
  m_tree->SetBranchAddress("filteredResiduals", &m_filteredResiduals);
  m_tree->SetBranchAddress("smoothedResiduals", &m_smoothedResiduals);

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

  /// Initial guess of the momentum at the IP
  m_tree->SetBranchAddress("originMomentumGuess", &m_originMomentumGuess);

  /// Initial guess of the vertex at the IP
  m_tree->SetBranchAddress("vertexGuess", &m_vertexGuess);

  /// KF predicted momentum at the IP
  m_tree->SetBranchAddress("originMomentumEst", &m_originMomentumEst);

  /// KF predicted vertex at the IP
  m_tree->SetBranchAddress("vertexEst", &m_vertexEst);

  // Chi2 and ndf of the fitted track
  m_tree->SetBranchAddress("chi2Predicted", &m_chi2Predicted);
  m_tree->SetBranchAddress("chi2Filtered", &m_chi2Filtered);
  m_tree->SetBranchAddress("chi2Smoothed", &m_chi2Smoothed);
  m_tree->SetBranchAddress("ndf", &m_ndf);

  // Track ID
  m_tree->SetBranchAddress("trackId", &m_trackId);

  // Event ID
  m_tree->SetBranchAddress("eventId", &m_eventId);

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
  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);

  m_outputSeedsGuess.initialize(m_cfg.outputSeedsGuess);
  m_outputTrackParametersGuess.initialize(m_cfg.outputTrackParametersGuess);

  m_outputSeedsEst.initialize(m_cfg.outputSeedsEst);
  m_outputTrackParametersEst.initialize(m_cfg.outputTrackParametersEst);
}

std::pair<std::size_t, std::size_t> RootTrackReader::availableEvents() const {
  return {std::get<0>(m_eventMap.front()), std::get<0>(m_eventMap.back()) + 1};
}

ProcessCode RootTrackReader::read(const AlgorithmContext& ctx) {
  auto it = std::ranges::find_if(m_eventMap, [&](const auto& a) {
    return std::get<0>(a) == ctx.eventNumber;
  });

  if (it == m_eventMap.end()) {
    // Explicitly warn if it happens for the first or last event as that might
    // indicate a human error
    if ((ctx.eventNumber == availableEvents().first) &&
        (ctx.eventNumber == availableEvents().second - 1)) {
      ACTS_WARNING("Reading empty event: " << ctx.eventNumber);
    } else {
      ACTS_DEBUG("Reading empty event: " << ctx.eventNumber);
    }

    m_outputSourceLinks(ctx, {});

    m_outputSeedsGuess(ctx, {});
    m_outputTrackParametersGuess(ctx, {});

    m_outputSeedsEst(ctx, {});
    m_outputTrackParametersEst(ctx, {});

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

  IndexSeeds seedsGuess{};
  std::vector<Acts::CurvilinearTrackParameters> trackParametersGuess{};

  IndexSeeds seedsEst{};
  std::vector<Acts::CurvilinearTrackParameters> trackParametersEst{};

  std::size_t eventId = std::get<0>(*it);
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
    if (m_originMomentumEst->P() < constraints.minAbsMomentumEst ||
        m_originMomentumEst->P() > constraints.maxAbsMomentumEst) {
      continue;
    }

    Acts::BoundMatrix ipCovGuess;
    Acts::BoundMatrix ipCovEst;
    for (std::size_t i = 0; i < Acts::eBoundSize; i++) {
      for (std::size_t j = 0; j < Acts::eBoundSize; j++) {
        ipCovGuess(i, j) = (*m_boundTrackCovGuess)(i, j);
        ipCovEst(i, j) = (*m_boundTrackCovEst)(i, j);
      }
    }

    std::vector<std::size_t> trackSourceLinkIndices{};
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
                                     geoId, eventId, sourceLinks.size());
      trackSourceLinkIndices.push_back(sourceLinks.size());
      sourceLinks.push_back(Acts::SourceLink(obsSourceLink));
    }

    // Initial guess
    Acts::ParticleHypothesis hypothesis =
        Acts::ParticleHypothesis(Acts::PdgParticle(m_pdgId));
    Acts::Vector4 vertexGuess(m_vertexGuess->X(), m_vertexGuess->Y(),
                              m_vertexGuess->Z(), 0);
    Acts::Vector3 ipDirectionGuess(m_originMomentumGuess->X(),
                                   m_originMomentumGuess->Y(),
                                   m_originMomentumGuess->Z());
    ipDirectionGuess.normalize();

    seedsGuess.emplace_back(trackSourceLinkIndices, trackParametersGuess.size(),
                            static_cast<int>(seedsGuess.size()));
    trackParametersGuess.emplace_back(vertexGuess, ipDirectionGuess,
                                      m_charge / m_originMomentumGuess->P(),
                                      ipCovGuess, hypothesis);

    // Estimated
    Acts::Vector4 vertexEst(m_vertexEst->X(), m_vertexEst->Y(),
                            m_vertexEst->Z(), 0);
    Acts::Vector3 ipDirectionEst(m_originMomentumEst->X(),
                                 m_originMomentumEst->Y(),
                                 m_originMomentumEst->Z());
    ipDirectionEst.normalize();

    seedsEst.emplace_back(trackSourceLinkIndices, trackParametersEst.size(),
                          static_cast<int>(seedsEst.size()));
    trackParametersEst.emplace_back(vertexEst, ipDirectionEst,
                                    m_charge / m_originMomentumEst->P(),
                                    ipCovEst, hypothesis);
  }

  ACTS_DEBUG("Read " << sourceLinks.size() << " source links");

  ACTS_DEBUG("Read " << seedsGuess.size() << " guess seeds");
  ACTS_DEBUG("Read " << trackParametersGuess.size()
                     << " guess track parameters");

  ACTS_DEBUG("Read " << seedsEst.size() << " estimated seeds");
  ACTS_DEBUG("Read " << trackParametersEst.size()
                     << " estimated track parameters");

  m_outputSourceLinks(ctx, std::move(sourceLinks));

  m_outputSeedsGuess(ctx, std::move(seedsGuess));
  m_outputTrackParametersGuess(ctx, std::move(trackParametersGuess));

  m_outputSeedsEst(ctx, std::move(seedsEst));
  m_outputTrackParametersEst(ctx, std::move(trackParametersEst));

  // Return success flag
  return ProcessCode::SUCCESS;
}
