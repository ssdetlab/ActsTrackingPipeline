#include "TrackingPipeline/Io/RootTrackWriter.hpp"

#include <Acts/Definitions/Algebra.hpp>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

RootTrackWriter::RootTrackWriter(const Config& config,
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

  // KF residuals with respect to the measurements
  m_tree->Branch("predictedResiduals", &m_predictedResiduals, buf_size,
                 split_lvl);
  m_tree->Branch("filteredResiduals", &m_filteredResiduals, buf_size,
                 split_lvl);
  m_tree->Branch("smoothedResiduals", &m_smoothedResiduals, buf_size,
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

  // Chi2 and ndf of the fitted track
  m_tree->Branch("chi2Predicted", &m_chi2Predicted, "chi2Predicted/D");
  m_tree->Branch("chi2Filtered", &m_chi2Filtered, "chi2Filtered/D");
  m_tree->Branch("chi2Smoothed", &m_chi2Smoothed, "chi2Smoothed/D");
  m_tree->Branch("ndf", &m_ndf, "ndf/I");

  // Track ID
  m_tree->Branch("trackId", &m_trackId, buf_size, split_lvl);

  // Event ID
  m_tree->Branch("eventId", &m_eventId, "eventId/I");

  // PDG ID
  m_tree->Branch("pdgId", &m_pdgId, "pdgId/I");

  // Charge
  m_tree->Branch("charge", &m_charge, "charge/I");

  //------------------------------------------------------------------
  // Initialize the data handles
  m_inputTracks.initialize(m_cfg.inputTracks);
}

ProcessCode RootTrackWriter::finalize() {
  if (m_file) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootTrackWriter::write(const AlgorithmContext& ctx) {
  auto inputTracks = m_inputTracks(ctx);

  std::lock_guard<std::mutex> lock(m_mutex);

  m_eventId = ctx.eventNumber;

  // Iterate over the fitted tracks
  for (std::size_t tid = 0; tid < inputTracks.tracks.size(); tid++) {
    // Get the track object and the track id
    const auto& track = inputTracks.tracks.getTrack(tid);

    // Get PDG id
    m_pdgId = inputTracks.ipParametersGuesses.at(tid)
                  .particleHypothesis()
                  .absolutePdg();

    // Get charge
    m_charge = inputTracks.ipParametersGuesses.at(tid).charge();

    // Guess momentum used for the KF fit
    m_ipMomentumGuess.SetPxPyPzE(
        inputTracks.ipParametersGuesses.at(tid).momentum().x(),
        inputTracks.ipParametersGuesses.at(tid).momentum().y(),
        inputTracks.ipParametersGuesses.at(tid).momentum().z(),
        std::hypot(inputTracks.ipParametersGuesses.at(tid).momentum().norm(),
                   inputTracks.ipParametersGuesses.at(tid)
                       .particleHypothesis()
                       .mass()));

    // Guess vertex used for the KF fit
    m_vertexGuess =
        TVector3(inputTracks.ipParametersGuesses.at(tid).position().x(),
                 inputTracks.ipParametersGuesses.at(tid).position().y(),
                 inputTracks.ipParametersGuesses.at(tid).position().z());

    // KF predicted momentum at the IP
    Acts::Vector3 pVec = track.momentum();
    double pMag = pVec.norm();
    m_ipMomentumEst.SetPxPyPzE(
        pVec.x(), pVec.y(), pVec.z(),
        std::hypot(pMag, inputTracks.ipParametersGuesses.at(tid)
                             .particleHypothesis()
                             .mass()));

    // KF predicted IP momentum error
    m_ipMomentumError =
        TVector3(std::sqrt(track.covariance().diagonal().head<4>()[2]),
                 std::sqrt(track.covariance().diagonal().head<4>()[3]), 0);

    // KF predicted vertex position
    Acts::Vector3 vertex = track.referenceSurface().localToGlobal(
        ctx.geoContext, Acts::Vector2{track.loc0(), track.loc1()},
        Acts::Vector3::UnitY());
    m_vertexEst = TVector3(vertex.x(), vertex.y(), vertex.z());

    // KF predicted vertex error
    Acts::Vector3 vertexError = {
        std::sqrt(track.covariance().diagonal().head<2>()[0]), 0,
        std::sqrt(track.covariance().diagonal().head<2>()[1])};
    m_vertexError = TVector3(vertexError.x(), vertexError.y(), vertexError.z());

    std::size_t nStates = track.nTrackStates();

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

    // KF residuals with respect to the measurements
    std::vector<TVector2> predictedResiduals;
    predictedResiduals.reserve(nStates);

    std::vector<TVector2> filteredResiduals;
    filteredResiduals.reserve(nStates);

    std::vector<TVector2> smoothedResiduals;
    smoothedResiduals.reserve(nStates);

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

    // Chi2 at different states
    double chi2Predicted = 0;
    double chi2Filtered = 0;
    double chi2Smoothed = 0;

    // Iterate over the track states
    for (const auto& state : track.trackStatesReversed()) {
      // Skip the states without meaningful information
      if (!state.hasProjector()) {
        continue;
      }
      if (!state.hasUncalibratedSourceLink()) {
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

      trackIds.push_back(inputTracks.trackIds.at(tid));

      // ---------------------------------------------
      // Track hit info

      // Get the measurements hit
      Acts::Vector2 hit = state.effectiveCalibrated();
      trackHitsLocal.push_back(TVector2(hit.x(), hit.y()));

      // Transform the hits to the global coordinates
      Acts::Vector3 hitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, hit, Acts::Vector3(1, 0, 0));

      // Covariance
      Acts::SquareMatrix2 measurementCov =
          state.effectiveCalibratedCovariance();

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

      // KF pulls with respect to the measurements
      Acts::Vector2 predictedPull =
          predictedDiag.cwiseProduct(hit - predictedHit);

      // Store the KF predicted hits
      predictedTrackHitsGlobal.push_back(TVector3(predictedHitGlobal.x(),
                                                  predictedHitGlobal.y(),
                                                  predictedHitGlobal.z()));

      // Store the residuals with respect to the measurements
      predictedResiduals.push_back(
          TVector2(predictedResidual.x(), predictedResidual.y()));

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

        Acts::Vector2 filteredResidual = hit - filteredHit;

        Acts::SquareMatrix2 filteredCov =
            state.effectiveProjector() * state.filteredCovariance() *
                state.effectiveProjector().transpose() -
            measurementCov;
        Acts::Vector2 filteredDiag =
            filteredCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

        Acts::Vector2 filteredPull =
            filteredDiag.cwiseProduct(hit - filteredHit);

        filteredTrackHitsGlobal.push_back(TVector3(filteredHitGlobal.x(),
                                                   filteredHitGlobal.y(),
                                                   filteredHitGlobal.z()));
        filteredResiduals.push_back(
            TVector2(filteredResidual.x(), filteredResidual.y()));

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

        Acts::Vector2 smoothedResidual = hit - smoothedHit;

        Acts::SquareMatrix2 smoothedCov =
            state.effectiveProjector() * state.smoothedCovariance() *
                state.effectiveProjector().transpose() -
            measurementCov;
        Acts::Vector2 smoothedDiag =
            smoothedCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

        Acts::Vector2 smoothedPull =
            smoothedDiag.cwiseProduct(hit - smoothedHit);

        smoothedTrackHitsGlobal.push_back(TVector3(smoothedHitGlobal.x(),
                                                   smoothedHitGlobal.y(),
                                                   smoothedHitGlobal.z()));
        smoothedResiduals.push_back(
            TVector2(smoothedResidual.x(), smoothedResidual.y()));

        smoothedPulls.push_back(TVector2(smoothedPull.x(), smoothedPull.y()));

        chi2Smoothed += smoothedPull.dot(smoothedPull);
      }
    }
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

    // KF residuals with respect to the measurements
    m_predictedResiduals = std::move(predictedResiduals);
    m_filteredResiduals = std::move(filteredResiduals);
    m_smoothedResiduals = std::move(smoothedResiduals);

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

    // Fill the tree
    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
