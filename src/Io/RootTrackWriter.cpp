#include "TrackingPipeline/Io/RootTrackWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"

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
  int bufSize = 32000;
  int splitLvl = 0;

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
  m_tree->Branch("smoothedTrackHitsGlobal", &m_smoothedTrackHitsGlobal, bufSize,
                 splitLvl);

  m_tree->Branch("predictedTrackHitsLocal", &m_predictedTrackHitsLocal, bufSize,
                 splitLvl);
  m_tree->Branch("filteredTrackHitsLocal", &m_filteredTrackHitsLocal, bufSize,
                 splitLvl);
  m_tree->Branch("smoothedTrackHitsLocal", &m_smoothedTrackHitsLocal, bufSize,
                 splitLvl);

  // KF residuals with respect to the measurements
  m_tree->Branch("predictedResiduals", &m_predictedResiduals, bufSize,
                 splitLvl);
  m_tree->Branch("filteredResiduals", &m_filteredResiduals, bufSize, splitLvl);
  m_tree->Branch("smoothedResiduals", &m_smoothedResiduals, bufSize, splitLvl);

  // KF pulls with respect to the measurements
  m_tree->Branch("predictedPulls", &m_predictedPulls, bufSize, splitLvl);
  m_tree->Branch("filteredPulls", &m_filteredPulls, bufSize, splitLvl);
  m_tree->Branch("smoothedPulls", &m_smoothedPulls, bufSize, splitLvl);

  /// Guessed bound track parameters
  m_tree->Branch("boundTrackParametersGuess", &m_boundTrackParametersGuess,
                 bufSize, splitLvl);
  m_tree->Branch("boundTrackCovGuess", &m_boundTrackCovGuess, bufSize,
                 splitLvl);

  /// KF predicted bound track parameters
  m_tree->Branch("boundTrackParametersEst", &m_boundTrackParametersEst, bufSize,
                 splitLvl);
  m_tree->Branch("boundTrackCovEst", &m_boundTrackCovEst, bufSize, splitLvl);

  /// Initial guess of the momentum at the IP
  m_tree->Branch("ipMomentumGuess", &m_ipMomentumGuess, bufSize, splitLvl);

  /// Initial guess of the vertex at the IP
  m_tree->Branch("vertexGuess", &m_vertexGuess, bufSize, splitLvl);

  /// KF predicted momentum at the IP
  m_tree->Branch("ipMomentumEst", &m_ipMomentumEst, bufSize, splitLvl);

  /// KF predicted vertex at the IP
  m_tree->Branch("vertexEst", &m_vertexEst, bufSize, splitLvl);

  // Chi2 and ndf of the fitted track
  m_tree->Branch("chi2Predicted", &m_chi2Predicted, bufSize, splitLvl);
  m_tree->Branch("chi2Filtered", &m_chi2Filtered, bufSize, splitLvl);
  m_tree->Branch("chi2Smoothed", &m_chi2Smoothed, bufSize, splitLvl);
  m_tree->Branch("ndf", &m_ndf, bufSize, splitLvl);

  // Chi2 per degree of freedom
  m_tree->Branch("chi2PerDoFPredicted", &m_chi2PerDoFPredicted, bufSize, splitLvl); //////////////
  m_tree->Branch("chi2PerDoFFiltered",  &m_chi2PerDoFFiltered,  bufSize, splitLvl); //////////////
  m_tree->Branch("chi2PerDoFSmoothed",  &m_chi2PerDoFSmoothed,  bufSize, splitLvl); //////////////

  // Track ID
  m_tree->Branch("trackId", &m_trackId, bufSize, splitLvl);

  // Event ID
  m_tree->Branch("eventId", &m_eventId, bufSize, splitLvl);

  // PDG ID
  m_tree->Branch("pdgId", &m_pdgId, bufSize, splitLvl);

  // Charge
  m_tree->Branch("charge", &m_charge, bufSize, splitLvl);

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
  const auto& inputTracks = m_inputTracks(ctx);

  std::lock_guard<std::mutex> lock(m_mutex);

  m_eventId = ctx.eventNumber;

  // Iterate over the fitted tracks
  for (std::size_t tid = 0; tid < inputTracks.tracks.size(); tid++) {
    // Get the track object and the track id
    const auto& track = inputTracks.tracks.getTrack(tid);

    m_trackId = inputTracks.trackIds.at(tid);
    std::size_t nStates = track.nTrackStates();

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

    m_smoothedTrackHitsGlobal.clear();
    m_smoothedTrackHitsGlobal.reserve(nStates);

    m_predictedTrackHitsLocal.clear();
    m_predictedTrackHitsLocal.reserve(nStates);

    m_filteredTrackHitsLocal.clear();
    m_filteredTrackHitsLocal.reserve(nStates);

    m_smoothedTrackHitsLocal.clear();
    m_smoothedTrackHitsLocal.reserve(nStates);

    // KF residuals with respect to the measurements
    m_predictedResiduals.clear();
    m_predictedResiduals.reserve(nStates);

    m_filteredResiduals.clear();
    m_filteredResiduals.reserve(nStates);

    m_smoothedResiduals.clear();
    m_smoothedResiduals.reserve(nStates);

    // KF pulls with respect to the measurements
    m_predictedPulls.clear();
    m_predictedPulls.reserve(nStates);

    m_filteredPulls.clear();
    m_filteredPulls.reserve(nStates);

    m_smoothedPulls.clear();
    m_smoothedPulls.reserve(nStates);

    // ----------------------------------------------
    // Guess track parameters

    const auto& ipParametersGuess = inputTracks.ipParametersGuesses.at(tid);

    // Guessed IP momentum
    const auto& ipMomentumGuess = ipParametersGuess.momentum();
    double particleMass = ipParametersGuess.particleHypothesis().mass();
    m_ipMomentumGuess.SetPxPyPzE(
        ipMomentumGuess.x(), ipMomentumGuess.y(), ipMomentumGuess.z(),
        std::hypot(ipMomentumGuess.norm(), particleMass));

    // Guessed vertex
    const auto& ipPositionGuess = ipParametersGuess.position(ctx.geoContext);
    m_vertexGuess =
        TVector3(ipPositionGuess.x(), ipPositionGuess.y(), ipPositionGuess.z());

    // Guessed bound track parameters
    Acts::BoundVector boundTrackParametersGuess =
        ipParametersGuess.parameters();

    TArrayD boundTrackParsGuessData(Acts::eBoundSize);
    for (std::size_t i = 0; i < Acts::eBoundSize; i++) {
      boundTrackParsGuessData[i] = boundTrackParametersGuess(i);
    }
    m_boundTrackParametersGuess.Use(0, Acts::eBoundSize,
                                    boundTrackParsGuessData.GetArray());

    // Guessed bound errors
    Acts::BoundMatrix boundTrackCovGuess =
        ipParametersGuess.covariance().value();
    TArrayD boundTrackCovGuessData(Acts::eBoundSize * Acts::eBoundSize);
    for (std::size_t i = 0; i < Acts::eBoundSize * Acts::eBoundSize; i++) {
      boundTrackCovGuessData[i] = boundTrackCovGuess(i);
    }
    m_boundTrackCovGuess.Use(Acts::eBoundSize, Acts::eBoundSize,
                             boundTrackCovGuessData.GetArray());

    // ----------------------------------------------
    // Estimated track parameters

    // Estimated IP momentum
    Acts::Vector3 ipMomentumEst = track.momentum();
    m_ipMomentumEst.SetPxPyPzE(ipMomentumEst.x(), ipMomentumEst.y(),
                               ipMomentumEst.z(),
                               std::hypot(ipMomentumEst.norm(), particleMass));

    // KF predicted vertex position
    Acts::Vector3 vertexEst = m_cfg.referenceSurface->localToGlobal(
        ctx.geoContext, {track.loc0(), track.loc1()},
        ipMomentumEst.normalized());
    m_vertexEst = TVector3(vertexEst.x(), vertexEst.y(), vertexEst.z());

    // KF predicted bound track parameters
    Acts::BoundVector boundTrackParametersEst = track.parameters();
    TArrayD boundTrackParsEstData(Acts::eBoundSize);
    for (std::size_t i = 0; i < Acts::eBoundSize; i++) {
      boundTrackParsEstData[i] = boundTrackParametersEst(i);
    }
    m_boundTrackParametersEst.Use(0, Acts::eBoundSize,
                                  boundTrackParsEstData.GetArray());

    // KF predicted bound errors
    Acts::BoundMatrix boundTrackCovEst = track.covariance();
    TArrayD boundTrackCovEstData(Acts::eBoundSize * Acts::eBoundSize);
    for (std::size_t i = 0; i < Acts::eBoundSize * Acts::eBoundSize; i++) {
      boundTrackCovEstData[i] = boundTrackCovEst(i);
    }
    m_boundTrackCovEst.Use(Acts::eBoundSize, Acts::eBoundSize,
                           boundTrackCovEstData.GetArray());

    // Get PDG id
    m_pdgId = ipParametersGuess.particleHypothesis().absolutePdg();

    // Get charge
    m_charge = ipParametersGuess.charge();

    // Get DoFs
    m_ndf = track.nDoF();

    // Iterate over the track states
    m_chi2Predicted = 0;
    m_chi2Filtered = 0;
    m_chi2Smoothed = 0;
    for (const auto& state : track.trackStatesReversed()) {
      // Skip the states without meaningful information
      if (!state.hasProjector()) {
        continue;
      }

      state.referenceSurface().geometryId();

      // Get the measurements source link
      auto ssl = state.getUncalibratedSourceLink().get<SimpleSourceLink>();

      m_geometryIds.push_back(ssl.geometryId().sensitive());
      TArrayD trackHitCovData(4);
      TMatrixD trackHitCov(2, 2);
      for (std::size_t i = 0; i < 4; i++) {
        trackHitCovData[i] = ssl.covariance()(i);
      }
      trackHitCov.Use(2, 2, trackHitCovData.GetArray());
      m_trackHitCovs.push_back(trackHitCov);

      // ---------------------------------------------
      // Track hit info

      // Get the measurements hit
      Acts::Vector2 hit = state.effectiveCalibrated();
      m_trackHitsLocal.emplace_back(hit.x(), hit.y());

      // Transform the hits to the global coordinates
      Acts::Vector3 hitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, hit, Acts::Vector3(1, 0, 0));

      // Covariance
      Acts::SquareMatrix2 measurementCov =
          state.effectiveCalibratedCovariance();

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
          ctx.geoContext, predictedHit, Acts::Vector3(1, 0, 0));

      // Get the residuals between the measurements and the predicted hits
      Acts::Vector2 predictedResidual = hit - predictedHit;

      // With respect to measurement
      Acts::SquareMatrix2 predictedCov = state.effectiveProjector() *
                                         state.predictedCovariance() *
                                         state.effectiveProjector().transpose();

      // Extract diagonals
      Acts::Vector2 predictedDiag =
          predictedCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

      // KF pulls with respect to the measurements
      Acts::Vector2 predictedPull =
          predictedDiag.cwiseProduct(hit - predictedHit);

      // Store the KF predicted hits
      m_predictedTrackHitsGlobal.emplace_back(predictedHitGlobal.x(),
                                              predictedHitGlobal.y(),
                                              predictedHitGlobal.z());

      // Store the residuals with respect to the measurements
      m_predictedResiduals.emplace_back(predictedResidual.x(),
                                        predictedResidual.y());

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

        m_filteredTrackHitsGlobal.emplace_back(filteredHitGlobal.x(),
                                               filteredHitGlobal.y(),
                                               filteredHitGlobal.z());

        m_filteredResiduals.emplace_back(filteredResidual.x(),
                                         filteredResidual.y());

        m_filteredPulls.emplace_back(filteredPull.x(), filteredPull.y());

        m_chi2Filtered += filteredPull.dot(filteredPull);
      }

      // ---------------------------------------------
      // Smoothed state info
      if (state.hasSmoothed()) {
        Acts::Vector2 smoothedHit =
            state.effectiveProjector() * state.smoothed();

        m_smoothedTrackHitsLocal.emplace_back(smoothedHit.x(), smoothedHit.y());
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

        m_smoothedTrackHitsGlobal.emplace_back(smoothedHitGlobal.x(),
                                               smoothedHitGlobal.y(),
                                               smoothedHitGlobal.z());

        m_smoothedResiduals.emplace_back(smoothedResidual.x(),
                                         smoothedResidual.y());

        m_smoothedPulls.emplace_back(smoothedPull.x(), smoothedPull.y());

        m_chi2Smoothed += smoothedPull.dot(smoothedPull);
      }
    }
    ///////////////////////////////////////////////////////////////////////////
    if (m_ndf > 0) {
      m_chi2PerDoFPredicted = m_chi2Predicted / static_cast<double>(m_ndf);
      m_chi2PerDoFFiltered  = m_chi2Filtered  / static_cast<double>(m_ndf);
      m_chi2PerDoFSmoothed  = m_chi2Smoothed  / static_cast<double>(m_ndf);
    } else {
      m_chi2PerDoFPredicted = 0.0;
      m_chi2PerDoFFiltered  = 0.0;
      m_chi2PerDoFSmoothed  = 0.0;
    }
    ///////////////////////////////////////////////////////////////////////////

    // Fill the tree
    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
