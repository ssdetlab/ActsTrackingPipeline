#include "TrackingPipeline/Io/RootSimTrackCandidateWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include <Acts/Utilities/Logger.hpp>

#include <algorithm>
#include <ranges>

#include <boost/mp11/algorithm.hpp>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

RootSimTrackCandidateWriter::RootSimTrackCandidateWriter(
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
  // Track tree branches
  int bufSize = 32000;
  int splitLvl = 0;

  // True hits
  m_tree->Branch("trueTrackHitsGlobal", &m_trueTrackHitsGlobal, bufSize,
                 splitLvl);
  m_tree->Branch("trueTrackHitsLocal", &m_trueTrackHitsLocal, bufSize,
                 splitLvl);
  m_tree->Branch("onSurfaceMomentum", &m_onSurfaceMomentum, bufSize, splitLvl);
  m_tree->Branch("isSignal", &m_isSignal, bufSize, splitLvl);

  // Measurement hits
  m_tree->Branch("trackHitsGlobal", &m_trackHitsGlobal, bufSize, splitLvl);
  m_tree->Branch("trackHitsLocal", &m_trackHitsLocal, bufSize, splitLvl);

  // Covariances of the track hits
  m_tree->Branch("trackHitsCovs", &m_trackHitCovs, bufSize, splitLvl);

  // Geometry ids of the track hits
  m_tree->Branch("geometryIds", &m_geometryIds, bufSize, splitLvl);

  // KF predicted track hits
  m_tree->Branch("predictedTrackHitsGlobal", &m_predictedTrackHitsGlobal,
                 bufSize, splitLvl);
  m_tree->Branch("filteredTrackHitsGlobal", &m_filteredTrackHitsGlobal, bufSize,
                 splitLvl);

  m_tree->Branch("predictedTrackHitsLocal", &m_predictedTrackHitsLocal, bufSize,
                 splitLvl);
  m_tree->Branch("filteredTrackHitsLocal", &m_filteredTrackHitsLocal, bufSize,
                 splitLvl);

  // KF residuals with respect to the true hits
  m_tree->Branch("truePredictedResiduals", &m_truePredictedResiduals, bufSize,
                 splitLvl);
  m_tree->Branch("trueFilteredResiduals", &m_trueFilteredResiduals, bufSize,
                 splitLvl);

  // KF residuals with respect to the measurements
  m_tree->Branch("predictedResiduals", &m_predictedResiduals, bufSize,
                 splitLvl);
  m_tree->Branch("filteredResiduals", &m_filteredResiduals, bufSize, splitLvl);

  // KF pulls with respect to the true hits
  m_tree->Branch("truePredictedPulls", &m_truePredictedPulls, bufSize,
                 splitLvl);
  m_tree->Branch("trueFilteredPulls", &m_trueFilteredPulls, bufSize, splitLvl);

  // KF pulls with respect to the measurements
  m_tree->Branch("predictedPulls", &m_predictedPulls, bufSize, splitLvl);
  m_tree->Branch("filteredPulls", &m_filteredPulls, bufSize, splitLvl);

  // Initial guess of the momentum at the IP
  m_tree->Branch("ipMomentumGuess", &m_ipMomentumGuess, bufSize, splitLvl);
  m_tree->Branch("vertexGuess", &m_vertexGuess, bufSize, splitLvl);

  // True momentum at the IP
  m_tree->Branch("ipMomentumTruth", &m_ipMomentumTruth, bufSize, splitLvl);
  m_tree->Branch("vertexTruth", &m_vertexTruth, bufSize, splitLvl);

  // Chi2 and ndf of the fitted track
  m_tree->Branch("chi2Predicted", &m_chi2Predicted, bufSize, splitLvl);
  m_tree->Branch("chi2Filtered", &m_chi2Filtered, bufSize, splitLvl);
  m_tree->Branch("ndf", &m_ndf, bufSize, splitLvl);

  // Track ID
  m_tree->Branch("stateTrackId", &m_stateTrackId, bufSize, splitLvl);
  m_tree->Branch("stateParentTrackId", &m_stateParentTrackId, bufSize,
                 splitLvl);
  m_tree->Branch("stateRunId", &m_stateRunId, bufSize, splitLvl);

  m_tree->Branch("trackId", &m_trackId, bufSize, splitLvl);
  m_tree->Branch("parentTrackId", &m_parentTrackId, bufSize, splitLvl);
  m_tree->Branch("runId", &m_runId, bufSize, splitLvl);

  // Event ID
  m_tree->Branch("eventId", &m_eventId, bufSize, splitLvl);

  // True track size
  m_tree->Branch("trueTrackSize", &m_trueTrackSize, bufSize, splitLvl);
  m_tree->Branch("trackInCandidateSize", &m_trackInCandidateSize, bufSize,
                 splitLvl);

  // PDG ID
  m_tree->Branch("pdgId", &m_pdgId, bufSize, splitLvl);

  // Charge
  m_tree->Branch("charge", &m_charge, bufSize, splitLvl);

  //------------------------------------------------------------------
  // Initialize the data handles
  m_inputTrackCandidates.initialize(m_cfg.inputTrackCandidates);
  m_inputTruthClusters.initialize(m_cfg.inputTruthClusters);
}

ProcessCode RootSimTrackCandidateWriter::finalize() {
  if (m_file) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootSimTrackCandidateWriter::write(const AlgorithmContext& ctx) {
  const auto& inputTrackCandidates = m_inputTrackCandidates(ctx);
  const auto& inputTruthClusters = m_inputTruthClusters(ctx);

  ACTS_DEBUG("Received " << inputTrackCandidates.size() << " track candidates");
  ACTS_DEBUG("Received " << inputTruthClusters.size() << " truth clusters");

  std::lock_guard<std::mutex> lock(m_mutex);

  // Collect true track statistics
  auto trueTrackIds =
      inputTruthClusters |
      std::views::filter([](const auto& cl) { return cl.isSignal; }) |
      std::views::transform([](const auto& cl) { return cl.truthHits; }) |
      std::views::join | std::views::transform([](const auto& hit) -> TrackID {
        return {hit.trackId, hit.parentTrackId, hit.runId};
      });

  m_eventId = ctx.eventNumber;

  // Iterate over the fitted tracks
  for (std::size_t tid = 0; tid < inputTrackCandidates.tracks.size(); tid++) {
    // Get the track object and the track id
    const auto& candidate = inputTrackCandidates.tracks.getTrack(tid);

    std::size_t nStates = candidate.nTrackStates();

    // Track hits from the true information
    m_trueTrackHitsGlobal.clear();
    m_trueTrackHitsGlobal.reserve(nStates);

    m_trueTrackHitsLocal.clear();
    m_trueTrackHitsLocal.reserve(nStates);

    m_onSurfaceMomentum.clear();
    m_onSurfaceMomentum.reserve(nStates);

    m_isSignal.clear();
    m_isSignal.reserve(nStates);

    // Covariances of the track hits
    m_trackHitCovs.clear();
    m_trackHitCovs.reserve(nStates);

    // Track hits gometry identifiers
    m_geometryIds.clear();
    m_geometryIds.reserve(nStates);

    // Track hits from the measurements
    m_trackHitsGlobal.clear();
    m_trackHitsGlobal.reserve(nStates);

    m_trackHitsLocal.clear();
    m_trackHitsLocal.reserve(nStates);

    // KF predicted track hits
    m_predictedTrackHitsGlobal.clear();
    m_predictedTrackHitsGlobal.reserve(nStates);

    m_filteredTrackHitsGlobal.clear();
    m_filteredTrackHitsGlobal.reserve(nStates);

    m_predictedTrackHitsLocal.clear();
    m_predictedTrackHitsLocal.reserve(nStates);

    m_filteredTrackHitsLocal.clear();
    m_filteredTrackHitsLocal.reserve(nStates);

    // KF residuals with respect to the true hits
    m_truePredictedResiduals.clear();
    m_truePredictedResiduals.reserve(nStates);

    m_trueFilteredResiduals.clear();
    m_trueFilteredResiduals.reserve(nStates);

    // KF residuals with respect to the measurements
    m_predictedResiduals.clear();
    m_predictedResiduals.reserve(nStates);

    m_filteredResiduals.clear();
    m_filteredResiduals.reserve(nStates);

    // KF pulls with respect to the true hits
    m_truePredictedPulls.clear();
    m_truePredictedPulls.reserve(nStates);

    m_trueFilteredPulls.clear();
    m_trueFilteredPulls.reserve(nStates);

    // KF pulls with respect to the measurements
    m_predictedPulls.clear();
    m_predictedPulls.reserve(nStates);

    m_filteredPulls.clear();
    m_filteredPulls.reserve(nStates);

    // Track ids
    m_stateTrackId.clear();
    m_stateTrackId.reserve(nStates);

    m_stateParentTrackId.clear();
    m_stateParentTrackId.reserve(nStates);

    m_stateRunId.clear();
    m_stateRunId.reserve(nStates);

    // Guess IP momentum
    const auto& ipParsGuess = inputTrackCandidates.ipParametersGuesses.at(tid);
    m_ipMomentumGuess.SetPxPyPzE(
        ipParsGuess.momentum().x(), ipParsGuess.momentum().y(),
        ipParsGuess.momentum().z(),
        std::hypot(ipParsGuess.momentum().norm(),
                   ipParsGuess.particleHypothesis().mass()));

    // Initial guess of the vertex
    m_vertexGuess.SetXYZ(ipParsGuess.position().x(), ipParsGuess.position().y(),
                         ipParsGuess.position().z());

    // Get PDG id
    m_pdgId = ipParsGuess.particleHypothesis().absolutePdg();

    // Get charge
    m_charge = ipParsGuess.charge();

    // Get DoFs
    m_ndf = candidate.nDoF();

    // Iterate over the track states
    std::map<TrackID, std::vector<int>> trackStateIds;
    for (const auto& state : candidate.trackStatesReversed()) {
      // Skip the states without meaningful information
      if (!state.hasProjector()) {
        continue;
      }

      // Get the measurements source link
      auto ssl = state.getUncalibratedSourceLink().get<SimpleSourceLink>();

      m_geometryIds.push_back(ssl.geometryId().sensitive());
      TArrayD data(4);
      TMatrixD cov(2, 2);
      for (std::size_t i = 0; i < 4; i++) {
        data[i] = ssl.covariance()(i);
      }
      cov.Use(2, 2, data.GetArray());
      m_trackHitCovs.push_back(cov);

      const auto& cluster = inputTruthClusters.at(ssl.index());
      m_isSignal.push_back(cluster.isSignal);

      // Get the true hit
      Acts::Vector2 trueHit;
      TrackID currentTrackId;
      TLorentzVector onSurfaceMom;
      if (cluster.truthHits.size() == 0 || !cluster.isSignal) {
        trueHit = ssl.parametersLoc();
        onSurfaceMom.SetPxPyPzE(0, 0, 0, 0);

        currentTrackId = std::make_tuple(-1, -1, -1);
      } else {
        auto sig = std::ranges::find_if(cluster.truthHits, [](const auto& hit) {
          return (hit.trackId == 1);
        });
        trueHit = sig->truthParameters.head<2>();

        double onSurfP =
            std::abs(1. / sig->truthParameters[Acts::eBoundQOverP]);
        onSurfaceMom.SetPxPyPzE(
            onSurfP * std::sin(sig->truthParameters[Acts::eBoundTheta]) *
                std::cos(sig->truthParameters[Acts::eBoundPhi]),
            onSurfP * std::sin(sig->truthParameters[Acts::eBoundTheta]) *
                std::sin(sig->truthParameters[Acts::eBoundPhi]),
            onSurfP * std::cos(sig->truthParameters[Acts::eBoundTheta]),
            std::hypot(onSurfP, sig->ipParameters.particleHypothesis().mass()));
        m_onSurfaceMomentum.push_back(onSurfaceMom);
        currentTrackId =
            std::make_tuple(sig->trackId, sig->parentTrackId, sig->runId);
      }

      m_stateTrackId.push_back(std::get<0>(currentTrackId));
      m_stateParentTrackId.push_back(std::get<1>(currentTrackId));
      m_stateRunId.push_back(std::get<2>(currentTrackId));

      trackStateIds[currentTrackId].push_back(ssl.index());
      m_trueTrackHitsLocal.emplace_back(trueHit.x(), trueHit.y());

      // ---------------------------------------------
      // Track hit info

      // Get the measurements hit
      const Acts::Vector2& hit = state.effectiveCalibrated();
      m_trackHitsLocal.emplace_back(hit.x(), hit.y());

      // Transform the hits to the global coordinates
      Acts::Vector3 trueHitGlobal =
          m_cfg.surfaceAccessor(Acts::SourceLink(cluster.sourceLink))
              ->localToGlobal(ctx.geoContext, trueHit, Acts::Vector3::UnitY());
      Acts::Vector3 hitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, hit, Acts::Vector3::UnitY());

      // Covariance
      Acts::SquareMatrix2 measurementCov =
          state.effectiveCalibratedCovariance();

      // Store the true hits
      m_trueTrackHitsGlobal.emplace_back(trueHitGlobal.x(), trueHitGlobal.y(),
                                         trueHitGlobal.z());

      // Store the measurements hits
      m_trackHitsGlobal.emplace_back(hitGlobal.x(), hitGlobal.y(),
                                     hitGlobal.z());

      // ---------------------------------------------
      // Predicted state info

      // Project onto the prediction space
      Acts::Vector2 predictedHit =
          state.effectiveProjector() * state.predicted();

      m_predictedTrackHitsLocal.emplace_back(predictedHit.x(),
                                             predictedHit.y());

      Acts::Vector3 predictedHitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, predictedHit, Acts::Vector3::UnitY());

      // Get the residuals between the true and the predicted hits
      Acts::Vector2 truePredictedResidual = trueHit - predictedHit;

      // Get the residuals between the measurements and the predicted hits
      Acts::Vector2 predictedResidual = hit - predictedHit;

      // With respect to truth
      Acts::SquareMatrix2 predictedCovTruth =
          measurementCov + state.effectiveProjector() *
                               state.predictedCovariance() *
                               state.effectiveProjector().transpose();

      // With respect to measurement
      Acts::SquareMatrix2 predictedCov = state.effectiveProjector() *
                                         state.predictedCovariance() *
                                         state.effectiveProjector().transpose();

      // Extract diagonals
      Acts::Vector2 predictedDiagTruth =
          predictedCovTruth.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

      Acts::Vector2 predictedDiag =
          predictedCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

      // KF pulls with respect to the true hits
      Acts::Vector2 truePredictedPull =
          predictedDiagTruth.cwiseProduct(trueHit - predictedHit);

      // KF pulls with respect to the measurements
      Acts::Vector2 predictedPull =
          predictedDiag.cwiseProduct(hit - predictedHit);

      // Store the KF predicted hits
      m_predictedTrackHitsGlobal.emplace_back(predictedHitGlobal.x(),
                                              predictedHitGlobal.y(),
                                              predictedHitGlobal.z());

      // Store the residuals with respect to the true hits
      m_truePredictedResiduals.emplace_back(truePredictedResidual.x(),
                                            truePredictedResidual.y());

      // Store the residuals with respect to the measurements
      m_predictedResiduals.emplace_back(predictedResidual.x(),
                                        predictedResidual.y());

      // Store the pulls with respect to the true hits
      m_truePredictedPulls.emplace_back(truePredictedPull.x(),
                                        truePredictedPull.y());

      // Store the pulls with respect to the measurements
      m_predictedPulls.emplace_back(predictedPull.x(), predictedPull.y());

      // Add to the track chi2
      m_chi2Predicted += predictedPull.dot(predictedPull);

      // ---------------------------------------------
      // Filtered state info
      if (state.hasFiltered()) {
        Acts::Vector2 filteredHit =
            state.effectiveProjector() * state.filtered();

        m_filteredTrackHitsLocal.emplace_back(filteredHit.x(), filteredHit.y());

        Acts::Vector3 filteredHitGlobal =
            state.referenceSurface().localToGlobal(ctx.geoContext, filteredHit,
                                                   Acts::Vector3::UnitY());

        Acts::Vector2 trueFilteredResidual = trueHit - filteredHit;
        Acts::Vector2 filteredResidual = hit - filteredHit;

        Acts::SquareMatrix2 filteredCovTruth =
            state.effectiveProjector() * state.filteredCovariance() *
            state.effectiveProjector().transpose();

        Acts::SquareMatrix2 filteredCov =
            state.effectiveProjector() * state.filteredCovariance() *
                state.effectiveProjector().transpose() -
            measurementCov;
        Acts::Vector2 filteredDiagTruth =
            filteredCovTruth.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
        Acts::Vector2 filteredDiag =
            filteredCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

        Acts::Vector2 trueFilteredPull =
            filteredDiagTruth.cwiseProduct(trueHit - filteredHit);
        Acts::Vector2 filteredPull =
            filteredDiag.cwiseProduct(hit - filteredHit);

        m_filteredTrackHitsGlobal.emplace_back(filteredHitGlobal.x(),
                                               filteredHitGlobal.y(),
                                               filteredHitGlobal.z());
        m_trueFilteredResiduals.emplace_back(trueFilteredResidual.x(),
                                             trueFilteredResidual.y());

        m_filteredResiduals.emplace_back(filteredResidual.x(),
                                         filteredResidual.y());

        m_trueFilteredPulls.emplace_back(trueFilteredPull.x(),
                                         trueFilteredPull.y());

        m_filteredPulls.emplace_back(filteredPull.x(), filteredPull.y());

        m_chi2Filtered += filteredPull.dot(filteredPull);
      }
    }

    // Matching degree is computed with respect
    // to the most often occuring signal track
    auto refTrackId = std::ranges::max_element(
        trackStateIds, [](const auto& pairA, const auto& pairB) {
          if (std::get<0>(pairA.first) == -1 &&
              std::get<0>(pairB.first) != -1) {
            return true;
          } else if (std::get<0>(pairA.first) != -1 &&
                     std::get<0>(pairB.first) == -1) {
            return false;
          } else {
            return pairA.second.size() < pairB.second.size();
          }
        });
    if (std::get<0>(refTrackId->first) == -1) {
      m_ipMomentumTruth.SetPxPyPzE(0, 0, 0, 0);
      m_vertexTruth.SetXYZ(0, 0, 0);

      m_trueTrackSize = 0;
      m_trackInCandidateSize = 0;

      m_trackId = 0;
      m_parentTrackId = 0;
      m_runId = 0;
    } else {
      // Get the true IP parameters
      int refIndex = refTrackId->second.at(0);
      const auto& cluster = inputTruthClusters.at(refIndex);
      auto pivotHit =
          std::ranges::find_if(cluster.truthHits, [&](const auto& hit) {
            TrackID id{hit.trackId, hit.parentTrackId, hit.runId};
            return (id == refTrackId->first);
          });

      const auto& trueIpPars = pivotHit->ipParameters;
      m_ipMomentumTruth.SetPxPyPzE(
          trueIpPars.momentum().x(), trueIpPars.momentum().y(),
          trueIpPars.momentum().z(),
          std::hypot(trueIpPars.absoluteMomentum(),
                     trueIpPars.particleHypothesis().mass()));

      m_vertexTruth.SetXYZ(trueIpPars.position().x(), trueIpPars.position().y(),
                           trueIpPars.position().z());

      // Compute matching degree
      m_trueTrackSize = std::ranges::count(trueTrackIds, refTrackId->first);
      m_trackInCandidateSize = refTrackId->second.size();

      // Store the dominant track id
      m_trackId = std::get<0>(refTrackId->first);
      m_parentTrackId = std::get<1>(refTrackId->first);
      m_runId = std::get<2>(refTrackId->first);
    }

    // Fill the tree
    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
