#include "TrackingPipeline/Io/RootSimTrackWriter.hpp"

#include <Acts/Definitions/Algebra.hpp>

#include <algorithm>
#include <ranges>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

RootSimTrackWriter::RootSimTrackWriter(const Config& config,
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
  // Track tree branches
  int buf_size = 32000;
  int split_lvl = 0;

  // True hits
  m_tree->Branch("trueTrackHitsGlobal", &m_trueTrackHitsGlobal, buf_size,
                 split_lvl);
  m_tree->Branch("trueTrackHitsLocal", &m_trueTrackHitsLocal, buf_size,
                 split_lvl);
  m_tree->Branch("onSurfaceMomentum", &m_onSurfaceMomentum, buf_size,
                 split_lvl);
  m_tree->Branch("isSignal", &m_isSignal, buf_size, split_lvl);

  // Measurement hits
  m_tree->Branch("trackHitsGlobal", &m_trackHitsGlobal, buf_size, split_lvl);
  m_tree->Branch("trackHitsLocal", &m_trackHitsLocal, buf_size, split_lvl);

  // Covariances of the track hits
  m_tree->Branch("trackHitsCovs", &m_trackHitCovs, buf_size, split_lvl);

  // Geometry ids of the track hits
  m_tree->Branch("geometryIds", &m_geometryIds, buf_size, split_lvl);

  // KF predicted track hits
  m_tree->Branch("predictedTrackHitsGlobal", &m_predictedTrackHitsGlobal,
                 buf_size, split_lvl);
  m_tree->Branch("filteredTrackHitsGlobal", &m_filteredTrackHitsGlobal,
                 buf_size, split_lvl);
  m_tree->Branch("smoothedTrackHitsGlobal", &m_smoothedTrackHitsGlobal,
                 buf_size, split_lvl);

  m_tree->Branch("predictedTrackHitsLocal", &m_predictedTrackHitsLocal,
                 buf_size, split_lvl);
  m_tree->Branch("filteredTrackHitsLocal", &m_filteredTrackHitsLocal, buf_size,
                 split_lvl);
  m_tree->Branch("smoothedTrackHitsLocal", &m_smoothedTrackHitsLocal, buf_size,
                 split_lvl);

  // KF residuals with respect to the true hits
  m_tree->Branch("truePredictedResiduals", &m_truePredictedResiduals, buf_size,
                 split_lvl);
  m_tree->Branch("trueFilteredResiduals", &m_trueFilteredResiduals, buf_size,
                 split_lvl);
  m_tree->Branch("trueSmoothedResiduals", &m_trueSmoothedResiduals, buf_size,
                 split_lvl);

  // KF residuals with respect to the measurements
  m_tree->Branch("predictedResiduals", &m_predictedResiduals, buf_size,
                 split_lvl);
  m_tree->Branch("filteredResiduals", &m_filteredResiduals, buf_size,
                 split_lvl);
  m_tree->Branch("smoothedResiduals", &m_smoothedResiduals, buf_size,
                 split_lvl);

  // KF pulls with respect to the true hits
  m_tree->Branch("truePredictedPulls", &m_truePredictedPulls, buf_size,
                 split_lvl);
  m_tree->Branch("trueFilteredPulls", &m_trueFilteredPulls, buf_size,
                 split_lvl);
  m_tree->Branch("trueSmoothedPulls", &m_trueSmoothedPulls, buf_size,
                 split_lvl);

  // KF pulls with respect to the measurements
  m_tree->Branch("predictedPulls", &m_predictedPulls, buf_size, split_lvl);
  m_tree->Branch("filteredPulls", &m_filteredPulls, buf_size, split_lvl);
  m_tree->Branch("smoothedPulls", &m_smoothedPulls, buf_size, split_lvl);

  // Initial guess of the momentum at the IP
  m_tree->Branch("ipMomentumGuess", &m_ipMomentumGuess);
  m_tree->Branch("vertexGuess", &m_vertexGuess);

  // KF predicted momentum at the IP
  m_tree->Branch("ipMomentumEst", &m_ipMomentumEst);
  m_tree->Branch("ipMomentumError", &m_ipMomentumError);
  m_tree->Branch("vertexEst", &m_vertexEst);
  m_tree->Branch("vertexError", &m_vertexError);

  // True momentum at the IP
  m_tree->Branch("ipMomentumTruth", &m_ipMomentumTruth);
  m_tree->Branch("vertexTruth", &m_vertexTruth);

  // Chi2 and ndf of the fitted track
  m_tree->Branch("chi2Predicted", &m_chi2Predicted, "chi2Predicted/D");
  m_tree->Branch("chi2Filtered", &m_chi2Filtered, "chi2Filtered/D");
  m_tree->Branch("chi2Smoothed", &m_chi2Smoothed, "chi2Smoothed/D");
  m_tree->Branch("ndf", &m_ndf, "ndf/I");

  // Matching degree between the true and the fitted track
  m_tree->Branch("matchingDegree", &m_matchingDegree, "matchingDegree/D");

  // Track ID
  m_tree->Branch("trackId", &m_trackId, buf_size, split_lvl);
  m_tree->Branch("parentTrackId", &m_parentTrackId, buf_size, split_lvl);
  m_tree->Branch("runId", &m_runId, buf_size, split_lvl);

  // Event ID
  m_tree->Branch("eventId", &m_eventId, "eventId/I");

  // True track size
  m_tree->Branch("trueTrackSize", &m_trueTrackSize, "trueTrackSize/I");

  // PDG ID
  m_tree->Branch("pdgId", &m_pdgId, "pdgId/I");

  // Charge
  m_tree->Branch("charge", &m_charge, "charge/I");

  //------------------------------------------------------------------
  // Initialize the data handles
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputTruthClusters.initialize(m_cfg.inputTruthClusters);
}

ProcessCode RootSimTrackWriter::finalize() {
  if (m_file) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootSimTrackWriter::write(const AlgorithmContext& ctx) {
  auto inputTracks = m_inputTracks(ctx);

  auto inputTruthClusters = m_inputTruthClusters(ctx);

  std::lock_guard<std::mutex> lock(m_mutex);

  // Collect true track statistics
  std::map<TrackID, int> trueTracksSig;

  // Collect true track statistics
  auto trueTrackIds =
      inputTruthClusters |
      std::views::filter([](const auto& cl) { return cl.isSignal; }) |
      std::views::transform([](const auto& cl) { return cl.truthHits; }) |
      std::views::join | std::views::transform([](const auto& hit) -> TrackID {
        return {hit.trackId, hit.parentTrackId, hit.runId};
      });

  for (const auto& id : trueTrackIds) {
    if (!trueTracksSig.contains(id)) {
      trueTracksSig[id] = 1;
    } else {
      trueTracksSig.at(id)++;
    }
  }

  m_eventId = ctx.eventNumber;

  // Iterate over the fitted tracks
  for (std::size_t tid = 0; tid < inputTracks.tracks.size(); tid++) {
    // Get the track object and the track id
    const auto& track = inputTracks.tracks.getTrack(tid);

    // Initial guess of the vertex
    m_vertexGuess =
        TVector3(inputTracks.ipParametersGuesses.at(tid).position().x(),
                 inputTracks.ipParametersGuesses.at(tid).position().y(),
                 inputTracks.ipParametersGuesses.at(tid).position().z());

    // KF predicted momentum at the IP
    Acts::Vector3 pVec = track.momentum();
    double pMag = pVec.norm();

    // KF predicted vertex position
    m_vertexEst = TVector3(track.loc0(),
                           track.referenceSurface().center(ctx.geoContext).y(),
                           -track.loc1());

    // KF predicted vertex error
    Acts::Vector3 vertexError = {
        std::sqrt(track.covariance().diagonal().head<2>()[0]), 0,
        std::sqrt(track.covariance().diagonal().head<2>()[1])};
    m_vertexError = TVector3(vertexError.x(), vertexError.y(), vertexError.z());

    std::size_t nStates = track.nTrackStates();

    // Track hits from the true information
    std::vector<TVector3> trueTrackHitsGlobal;
    trueTrackHitsGlobal.reserve(nStates);

    std::vector<TVector2> trueTrackHitsLocal;
    trueTrackHitsLocal.reserve(nStates);

    std::vector<TLorentzVector> onSurfaceMomenta;
    onSurfaceMomenta.reserve(nStates);

    std::vector<int> isSignal;
    isSignal.reserve(nStates);

    // Covariances of the track hits
    std::vector<TMatrixD> trackHitCovs;
    trackHitCovs.reserve(nStates);

    // Track hits gometry identifiers
    std::vector<int> geometryIds;
    geometryIds.reserve(nStates);

    // Track hits from the measurements
    std::vector<TVector3> trackHitsGlobal;
    trackHitsGlobal.reserve(nStates);

    std::vector<TVector2> trackHitsLocal;
    trackHitsLocal.reserve(nStates);

    // KF predicted track hits
    std::vector<TVector3> predictedTrackHitsGlobal;
    predictedTrackHitsGlobal.reserve(nStates);

    std::vector<TVector3> filteredTrackHitsGlobal;
    filteredTrackHitsGlobal.reserve(nStates);

    std::vector<TVector3> smoothedTrackHitsGlobal;
    smoothedTrackHitsGlobal.reserve(nStates);

    std::vector<TVector2> predictedTrackHitsLocal;
    predictedTrackHitsLocal.reserve(nStates);

    std::vector<TVector2> filteredTrackHitsLocal;
    filteredTrackHitsLocal.reserve(nStates);

    std::vector<TVector2> smoothedTrackHitsLocal;
    smoothedTrackHitsLocal.reserve(nStates);

    // KF residuals with respect to the true hits
    std::vector<TVector2> truePredictedResiduals;
    truePredictedResiduals.reserve(nStates);

    std::vector<TVector2> trueFilteredResiduals;
    trueFilteredResiduals.reserve(nStates);

    std::vector<TVector2> trueSmoothedResiduals;
    trueSmoothedResiduals.reserve(nStates);

    // KF residuals with respect to the measurements
    std::vector<TVector2> predictedResiduals;
    predictedResiduals.reserve(nStates);

    std::vector<TVector2> filteredResiduals;
    filteredResiduals.reserve(nStates);

    std::vector<TVector2> smoothedResiduals;
    smoothedResiduals.reserve(nStates);

    // KF pulls with respect to the true hits
    std::vector<TVector2> truePredictedPulls;
    truePredictedPulls.reserve(nStates);

    std::vector<TVector2> trueFilteredPulls;
    trueFilteredPulls.reserve(nStates);

    std::vector<TVector2> trueSmoothedPulls;
    trueSmoothedPulls.reserve(nStates);

    // KF pulls with respect to the measurements
    std::vector<TVector2> predictedPulls;
    predictedPulls.reserve(nStates);

    std::vector<TVector2> filteredPulls;
    filteredPulls.reserve(nStates);

    std::vector<TVector2> smoothedPulls;
    smoothedPulls.reserve(nStates);

    // Track ids
    std::vector<int> trackIds;
    trackIds.reserve(nStates);

    std::vector<int> parentTrackIds;
    parentTrackIds.reserve(nStates);

    std::vector<int> runIds;
    runIds.reserve(nStates);

    // Flag indicating how many hits are matched
    // between the true and the fitted track
    double matchingDegree = 0;

    // Number of measurements in a true track
    int trueTrackSize = 0;

    // True particle PDG id
    int pdgId = 0;

    // True particle charge
    int charge = 0;

    // Chi2 at different states
    double chi2Predicted = 0;
    double chi2Filtered = 0;
    double chi2Smoothed = 0;

    // Iterate over the track states
    std::map<TrackID, std::vector<int>> trackStateIds;
    for (const auto& state : track.trackStatesReversed()) {
      // Skip the states without meaningful information
      if (!state.hasProjector()) {
        continue;
      }

      // Get the measurements source link
      auto sl = state.getUncalibratedSourceLink();
      auto ssl = sl.get<SimpleSourceLink>();

      geometryIds.push_back(ssl.geometryId().sensitive());
      TArrayD data(4);
      TMatrixD cov(2, 2);
      for (std::size_t i = 0; i < 4; i++) {
        data[i] = ssl.covariance()(i);
      }
      cov.Use(2, 2, data.GetArray());
      trackHitCovs.push_back(cov);

      auto cluster = inputTruthClusters.at(ssl.index());
      isSignal.push_back(cluster.isSignal);

      // Get the true hit
      Acts::Vector2 trueHit;
      TrackID currentTrackId;
      TLorentzVector onSurfaceMom;
      if (cluster.truthHits.size() == 0 || !cluster.isSignal) {
        trueHit = ssl.parameters();
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
        onSurfaceMomenta.push_back(onSurfaceMom);
        currentTrackId =
            std::make_tuple(sig->trackId, sig->parentTrackId, sig->runId);
      }

      trackIds.push_back(std::get<0>(currentTrackId));
      parentTrackIds.push_back(std::get<1>(currentTrackId));
      runIds.push_back(std::get<2>(currentTrackId));

      if (!trackStateIds.contains(currentTrackId)) {
        trackStateIds[currentTrackId] = {ssl.index()};
      } else {
        trackStateIds.at(currentTrackId).push_back(ssl.index());
      }
      trueTrackHitsLocal.push_back(TVector2(trueHit.x(), trueHit.y()));

      // ---------------------------------------------
      // Track hit info

      // Get the true source link
      auto trueSl = Acts::SourceLink(cluster.sourceLink);

      // Get the measurements hit
      Acts::Vector2 hit = state.effectiveCalibrated();
      trackHitsLocal.push_back(TVector2(hit.x(), hit.y()));

      // Transform the hits to the global coordinates
      Acts::Vector3 trueHitGlobal =
          m_cfg.surfaceAccessor(trueSl)->localToGlobal(ctx.geoContext, trueHit,
                                                       Acts::Vector3(1, 0, 0));
      Acts::Vector3 hitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, hit, Acts::Vector3(1, 0, 0));

      // Covariance
      Acts::SquareMatrix2 measurementCov =
          state.effectiveCalibratedCovariance();

      // Store the true hits
      trueTrackHitsGlobal.push_back(
          TVector3(trueHitGlobal.x(), trueHitGlobal.y(), trueHitGlobal.z()));

      // Store the measurements hits
      trackHitsGlobal.push_back(
          TVector3(hitGlobal.x(), hitGlobal.y(), hitGlobal.z()));

      // ---------------------------------------------
      // Predicted state info

      // Project onto the prediction space
      Acts::Vector2 predictedHit =
          state.effectiveProjector() * state.predicted();

      predictedTrackHitsLocal.push_back(
          TVector2(predictedHit.x(), predictedHit.y()));

      Acts::Vector3 predictedHitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, predictedHit, Acts::Vector3(1, 0, 0));

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
      predictedTrackHitsGlobal.push_back(TVector3(predictedHitGlobal.x(),
                                                  predictedHitGlobal.y(),
                                                  predictedHitGlobal.z()));

      // Store the residuals with respect to the true hits
      truePredictedResiduals.push_back(
          TVector2(truePredictedResidual.x(), truePredictedResidual.y()));

      // Store the residuals with respect to the measurements
      predictedResiduals.push_back(
          TVector2(predictedResidual.x(), predictedResidual.y()));

      // Store the pulls with respect to the true hits
      truePredictedPulls.push_back(
          TVector2(truePredictedPull.x(), truePredictedPull.y()));

      // Store the pulls with respect to the measurements
      predictedPulls.push_back(TVector2(predictedPull.x(), predictedPull.y()));

      // Add to the track chi2
      chi2Predicted += predictedPull.dot(predictedPull);

      // ---------------------------------------------
      // Filtered state info
      if (state.hasFiltered()) {
        Acts::Vector2 filteredHit =
            state.effectiveProjector() * state.filtered();

        filteredTrackHitsLocal.push_back(
            TVector2(filteredHit.x(), filteredHit.y()));

        Acts::Vector3 filteredHitGlobal =
            state.referenceSurface().localToGlobal(ctx.geoContext, filteredHit,
                                                   Acts::Vector3(1, 0, 0));

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

        filteredTrackHitsGlobal.push_back(TVector3(filteredHitGlobal.x(),
                                                   filteredHitGlobal.y(),
                                                   filteredHitGlobal.z()));
        trueFilteredResiduals.push_back(
            TVector2(trueFilteredResidual.x(), trueFilteredResidual.y()));

        filteredResiduals.push_back(
            TVector2(filteredResidual.x(), filteredResidual.y()));

        trueFilteredPulls.push_back(
            TVector2(trueFilteredPull.x(), trueFilteredPull.y()));

        filteredPulls.push_back(TVector2(filteredPull.x(), filteredPull.y()));

        chi2Filtered += filteredPull.dot(filteredPull);
      }

      // ---------------------------------------------
      // Smoothed state info
      if (state.hasSmoothed()) {
        Acts::Vector2 smoothedHit =
            state.effectiveProjector() * state.smoothed();

        smoothedTrackHitsLocal.push_back(
            TVector2(smoothedHit.x(), smoothedHit.y()));
        Acts::Vector3 smoothedHitGlobal =
            state.referenceSurface().localToGlobal(ctx.geoContext, smoothedHit,
                                                   Acts::Vector3(1, 0, 0));

        Acts::Vector2 trueSmoothedResidual = trueHit - smoothedHit;
        Acts::Vector2 smoothedResidual = hit - smoothedHit;

        Acts::SquareMatrix2 smoothedCovTruth =
            state.effectiveProjector() * state.smoothedCovariance() *
            state.effectiveProjector().transpose();

        Acts::SquareMatrix2 smoothedCov =
            state.effectiveProjector() * state.smoothedCovariance() *
                state.effectiveProjector().transpose() -
            measurementCov;
        Acts::Vector2 smoothedDiagTruth =
            smoothedCovTruth.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
        Acts::Vector2 smoothedDiag =
            smoothedCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

        Acts::Vector2 trueSmoothedPull =
            smoothedDiagTruth.cwiseProduct(trueHit - smoothedHit);
        Acts::Vector2 smoothedPull =
            smoothedDiag.cwiseProduct(hit - smoothedHit);

        smoothedTrackHitsGlobal.push_back(TVector3(smoothedHitGlobal.x(),
                                                   smoothedHitGlobal.y(),
                                                   smoothedHitGlobal.z()));
        trueSmoothedResiduals.push_back(
            TVector2(trueSmoothedResidual.x(), trueSmoothedResidual.y()));

        smoothedResiduals.push_back(
            TVector2(smoothedResidual.x(), smoothedResidual.y()));

        trueSmoothedPulls.push_back(
            TVector2(trueSmoothedPull.x(), trueSmoothedPull.y()));

        smoothedPulls.push_back(TVector2(smoothedPull.x(), smoothedPull.y()));

        chi2Smoothed += smoothedPull.dot(smoothedPull);
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
      matchingDegree = 0;
      trueTrackSize = 0;
    } else {
      // Get the true IP parameters
      int refIndex = refTrackId->second.at(0);
      auto cluster = inputTruthClusters.at(refIndex);
      auto pivotHit =
          std::ranges::find_if(cluster.truthHits, [&](const auto& hit) {
            TrackID id{hit.trackId, hit.parentTrackId, hit.runId};
            return (id == refTrackId->first);
          });

      m_ipMomentumTruth.SetPxPyPzE(
          pivotHit->ipParameters.momentum().x(),
          pivotHit->ipParameters.momentum().y(),
          pivotHit->ipParameters.momentum().z(),
          std::hypot(pivotHit->ipParameters.absoluteMomentum(),
                     pivotHit->ipParameters.particleHypothesis().mass()));

      // Compute matching degree
      trueTrackSize = std::ranges::count(trueTrackIds, refTrackId->first);

      matchingDegree =
          refTrackId->second.size() / static_cast<double>(trueTrackSize);

      // Get PDG id
      pdgId = pivotHit->ipParameters.particleHypothesis().absolutePdg();

      // Get charge
      charge = pivotHit->ipParameters.charge();

      // Guess IP momentum
      m_ipMomentumGuess.SetPxPyPzE(
          inputTracks.ipParametersGuesses.at(tid).momentum().x(),
          inputTracks.ipParametersGuesses.at(tid).momentum().y(),
          inputTracks.ipParametersGuesses.at(tid).momentum().z(),
          std::hypot(inputTracks.ipParametersGuesses.at(tid).momentum().norm(),
                     pivotHit->ipParameters.particleHypothesis().mass()));

      // Estimated IP momentum
      m_ipMomentumEst.SetPxPyPzE(
          pVec.x(), pVec.y(), pVec.z(),
          std::hypot(pMag, pivotHit->ipParameters.particleHypothesis().mass()));

      m_charge = pivotHit->ipParameters.charge();

      m_pdgId = pivotHit->ipParameters.particleHypothesis().absolutePdg();

      // KF predicted IP momentum error
      m_ipMomentumError =
          TVector3(std::sqrt(track.covariance().diagonal().head<4>()[2]),
                   std::sqrt(track.covariance().diagonal().head<4>()[3]), 0);
    }

    // True hits
    m_trueTrackHitsGlobal = std::move(trueTrackHitsGlobal);
    m_trueTrackHitsLocal = std::move(trueTrackHitsLocal);
    m_onSurfaceMomentum = std::move(onSurfaceMomenta);
    m_isSignal = std::move(isSignal);

    // Measurement hits
    m_trackHitsGlobal = std::move(trackHitsGlobal);
    m_trackHitsLocal = std::move(trackHitsLocal);

    // Covariances of the track hits
    m_trackHitCovs = std::move(trackHitCovs);

    // Geometry ids of the track hits
    m_geometryIds = std::move(geometryIds);

    // KF predicted track hits
    m_predictedTrackHitsGlobal = std::move(predictedTrackHitsGlobal);
    m_filteredTrackHitsGlobal = std::move(filteredTrackHitsGlobal);
    m_smoothedTrackHitsGlobal = std::move(smoothedTrackHitsGlobal);

    m_predictedTrackHitsLocal = std::move(predictedTrackHitsLocal);
    m_filteredTrackHitsLocal = std::move(filteredTrackHitsLocal);
    m_smoothedTrackHitsLocal = std::move(smoothedTrackHitsLocal);

    // KF residuals with respect to the true hits
    m_truePredictedResiduals = std::move(truePredictedResiduals);
    m_trueFilteredResiduals = std::move(trueFilteredResiduals);
    m_trueSmoothedResiduals = std::move(trueSmoothedResiduals);

    // KF residuals with respect to the measurements
    m_predictedResiduals = std::move(predictedResiduals);
    m_filteredResiduals = std::move(filteredResiduals);
    m_smoothedResiduals = std::move(smoothedResiduals);

    // KF pulls with respect to the true hits
    m_truePredictedPulls = std::move(truePredictedPulls);
    m_trueFilteredPulls = std::move(trueFilteredPulls);
    m_trueSmoothedPulls = std::move(trueSmoothedPulls);

    // KF pulls with respect to the measurements
    m_predictedPulls = std::move(predictedPulls);
    m_filteredPulls = std::move(filteredPulls);
    m_smoothedPulls = std::move(smoothedPulls);

    // Chi2 of the tracks from different
    // KF states
    m_chi2Predicted = chi2Predicted;
    m_chi2Filtered = chi2Filtered;
    m_chi2Smoothed = chi2Smoothed;

    // Number of degrees of freedom
    m_ndf = track.nDoF();

    // Track Id
    m_trackId = std::move(trackIds);
    m_parentTrackId = std::move(parentTrackIds);
    m_runId = std::move(runIds);

    // Matching degree
    m_matchingDegree = matchingDegree;

    // True track size
    m_trueTrackSize = trueTrackSize;

    // Fill the tree
    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
