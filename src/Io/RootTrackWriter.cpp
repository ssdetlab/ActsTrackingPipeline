#include "TrackingPipeline/Io/RootTrackWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"

#include <cstddef>
#include <stdexcept>

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

  // Measurement hits in the surface frame
  m_tree->Branch("trackHitsLocal", &m_trackHitsLocal, bufSize, splitLvl);

  // Measurement hits in the global frame
  m_tree->Branch("trackHitsGlobal", &m_trackHitsGlobal, bufSize, splitLvl);

  // Direction measurements in the track frame
  m_tree->Branch("onSurfaceTrackDirection", &m_onSurfaceTrackDirection, bufSize,
                 splitLvl);

  // Covariances of the track hits
  m_tree->Branch("trackHitCovs", &m_trackHitCovs, bufSize, splitLvl);

  // Covariances of the track agnles
  m_tree->Branch("trackAngleCovs", &m_trackAngleCovs, bufSize, splitLvl);

  // Geometry ids of the track hits
  m_tree->Branch("geometryIds", &m_geometryIds, bufSize, splitLvl);

  // KF predicted track hits in the surface frame
  m_tree->Branch("predictedTrackHitsLocal", &m_predictedTrackHitsLocal, bufSize,
                 splitLvl);
  m_tree->Branch("filteredTrackHitsLocal", &m_filteredTrackHitsLocal, bufSize,
                 splitLvl);
  m_tree->Branch("smoothedTrackHitsLocal", &m_smoothedTrackHitsLocal, bufSize,
                 splitLvl);

  // KF predicted track hits in the global frame
  m_tree->Branch("predictedTrackHitsGlobal", &m_predictedTrackHitsGlobal,
                 bufSize, splitLvl);
  m_tree->Branch("filteredTrackHitsGlobal", &m_filteredTrackHitsGlobal, bufSize,
                 splitLvl);
  m_tree->Branch("smoothedTrackHitsGlobal", &m_smoothedTrackHitsGlobal, bufSize,
                 splitLvl);

  // KF predicted on surface momenta in the track frame
  m_tree->Branch("predictedOnSurfaceMomentum", &m_predictedOnSurfaceMomentum,
                 bufSize, splitLvl);
  m_tree->Branch("filteredOnSurfaceMomentum", &m_filteredOnSurfaceMomentum,
                 bufSize, splitLvl);
  m_tree->Branch("smoothedOnSurfaceMomentum", &m_smoothedOnSurfaceMomentum,
                 bufSize, splitLvl);

  // KF residuals with respect to the measurements
  m_tree->Branch("predictedHitResiduals", &m_predictedHitResiduals, bufSize,
                 splitLvl);
  m_tree->Branch("filteredHitResiduals", &m_filteredHitResiduals, bufSize,
                 splitLvl);
  m_tree->Branch("smoothedHitResiduals", &m_smoothedHitResiduals, bufSize,
                 splitLvl);

  m_tree->Branch("predictedAngleResiduals", &m_predictedAngleResiduals, bufSize,
                 splitLvl);
  m_tree->Branch("filteredAngleResiduals", &m_filteredAngleResiduals, bufSize,
                 splitLvl);
  m_tree->Branch("smoothedAngleResiduals", &m_smoothedAngleResiduals, bufSize,
                 splitLvl);

  // KF pulls with respect to the measurements
  m_tree->Branch("predictedHitPulls", &m_predictedHitPulls, bufSize, splitLvl);
  m_tree->Branch("filteredHitPulls", &m_filteredHitPulls, bufSize, splitLvl);
  m_tree->Branch("smoothedHitPulls", &m_smoothedHitPulls, bufSize, splitLvl);

  m_tree->Branch("predictedAnglePulls", &m_predictedAnglePulls, bufSize,
                 splitLvl);
  m_tree->Branch("filteredAnglePulls", &m_filteredAnglePulls, bufSize,
                 splitLvl);
  m_tree->Branch("smoothedAnglePulls", &m_smoothedAnglePulls, bufSize,
                 splitLvl);

  // Guessed bound track parameters
  m_tree->Branch("boundTrackParametersGuess", &m_boundTrackParametersGuess,
                 bufSize, splitLvl);
  m_tree->Branch("boundTrackCovGuess", &m_boundTrackCovGuess, bufSize,
                 splitLvl);

  // KF predicted bound track parameters
  m_tree->Branch("boundTrackParametersEst", &m_boundTrackParametersEst, bufSize,
                 splitLvl);
  m_tree->Branch("boundTrackCovEst", &m_boundTrackCovEst, bufSize, splitLvl);

  // Initial guess of the momentum at the IP
  m_tree->Branch("originMomentumGuess", &m_originMomentumGuess, bufSize,
                 splitLvl);

  // Initial guess of the vertex at the IP
  m_tree->Branch("vertexGuess", &m_vertexGuess, bufSize, splitLvl);

  // KF predicted momentum at the IP
  m_tree->Branch("originMomentumEst", &m_originMomentumEst, bufSize, splitLvl);

  // KF predicted vertex at the IP
  m_tree->Branch("vertexEst", &m_vertexEst, bufSize, splitLvl);

  // Chi2 of the track with respect ot the measurement
  m_tree->Branch("chi2Predicted", &m_chi2Predicted, bufSize, splitLvl);
  m_tree->Branch("chi2Filtered", &m_chi2Filtered, bufSize, splitLvl);
  m_tree->Branch("chi2Smoothed", &m_chi2Smoothed, bufSize, splitLvl);

  // Number of degrees of freedom of the track
  m_tree->Branch("ndf", &m_ndf, bufSize, splitLvl);

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
  m_inputTrackContainer.initialize(m_cfg.inputTrackContainer);
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputTrackParametersGuesses.initialize(m_cfg.inputTrackParametersGuesses);
}

ProcessCode RootTrackWriter::finalize() {
  if (m_file != nullptr) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootTrackWriter::write(const AlgorithmContext& ctx) {
  const auto& inputTrackContainer = m_inputTrackContainer(ctx);
  const auto& inputTracks = m_inputTracks(ctx);
  const auto& inputTrackParametersGuesses = m_inputTrackParametersGuesses(ctx);

  ACTS_DEBUG("Received " << inputTrackContainer.size() << " input tracks");
  ACTS_DEBUG("Received " << inputTracks.size() << " input track indices");
  ACTS_DEBUG("Received " << inputTrackParametersGuesses.size()
                         << " input track parameters guesses");

  std::lock_guard<std::mutex> lock(m_mutex);

  // Event ID
  m_eventId = ctx.eventNumber;

  // Iterate over the fitted tracks
  for (std::size_t idx = 0; idx < inputTracks.size(); idx++) {
    // Get track indices
    const auto& trackIndices = inputTracks.at(idx);

    // Get the track object and the track id
    const auto& track = inputTrackContainer.getTrack(trackIndices.trackIndex);

    std::size_t nStates = track.nTrackStates();

    // Measurement hits in the surface frame
    m_trackHitsLocal.clear();
    m_trackHitsLocal.reserve(nStates);

    // Measurement hits in the global frame
    m_trackHitsGlobal.clear();
    m_trackHitsGlobal.reserve(nStates);

    // Direction measurements in the track frame
    m_onSurfaceTrackDirection.clear();
    m_onSurfaceTrackDirection.reserve(nStates);

    // Covariances of the track hits
    m_trackHitCovs.clear();
    m_trackHitCovs.reserve(nStates);

    // Covariances of the track agnles
    m_trackAngleCovs.clear();
    m_trackAngleCovs.reserve(nStates);

    // Geometry ids of the track hits
    m_geometryIds.clear();
    m_geometryIds.reserve(nStates);

    // KF predicted track hits in the surface frame
    m_predictedTrackHitsLocal.clear();
    m_predictedTrackHitsLocal.reserve(nStates);

    m_filteredTrackHitsLocal.clear();
    m_filteredTrackHitsLocal.reserve(nStates);

    m_smoothedTrackHitsLocal.clear();
    m_smoothedTrackHitsLocal.reserve(nStates);

    // KF predicted track hits in the global frame
    m_predictedTrackHitsGlobal.clear();
    m_predictedTrackHitsGlobal.reserve(nStates);

    m_filteredTrackHitsGlobal.clear();
    m_filteredTrackHitsGlobal.reserve(nStates);

    m_smoothedTrackHitsGlobal.clear();
    m_smoothedTrackHitsGlobal.reserve(nStates);

    // KF predicted on surface momenta in the track frame
    m_predictedOnSurfaceMomentum.clear();
    m_predictedOnSurfaceMomentum.reserve(nStates);

    m_filteredOnSurfaceMomentum.clear();
    m_filteredOnSurfaceMomentum.reserve(nStates);

    m_smoothedOnSurfaceMomentum.clear();
    m_smoothedOnSurfaceMomentum.reserve(nStates);

    // KF residuals with respect to the measurements
    m_predictedHitResiduals.clear();
    m_predictedHitResiduals.reserve(nStates);

    m_filteredHitResiduals.clear();
    m_filteredHitResiduals.reserve(nStates);

    m_smoothedHitResiduals.clear();
    m_smoothedHitResiduals.reserve(nStates);

    m_predictedAngleResiduals.clear();
    m_predictedAngleResiduals.reserve(nStates);

    m_filteredAngleResiduals.clear();
    m_filteredAngleResiduals.reserve(nStates);

    m_smoothedAngleResiduals.clear();
    m_smoothedAngleResiduals.reserve(nStates);

    // KF pulls with respect to the measurements
    m_predictedHitPulls.clear();
    m_predictedHitPulls.reserve(nStates);

    m_filteredHitPulls.clear();
    m_filteredHitPulls.reserve(nStates);

    m_smoothedHitPulls.clear();
    m_smoothedHitPulls.reserve(nStates);

    m_predictedAnglePulls.clear();
    m_predictedAnglePulls.reserve(nStates);

    m_filteredAnglePulls.clear();
    m_filteredAnglePulls.reserve(nStates);

    m_smoothedAnglePulls.clear();
    m_smoothedAnglePulls.reserve(nStates);

    // ----------------------------------------------
    // Guess track parameters

    const auto& originParametersGuess =
        inputTrackParametersGuesses.at(trackIndices.originParametersGuessIndex);

    // Guessed origin momentum
    const auto& originMomentumGuess = originParametersGuess.momentum();
    double particleMass = originParametersGuess.particleHypothesis().mass();
    m_originMomentumGuess.SetPxPyPzE(
        originMomentumGuess.x(), originMomentumGuess.y(),
        originMomentumGuess.z(),
        std::hypot(originMomentumGuess.norm(), particleMass));

    // Guessed vertex
    const auto& originPositionGuess =
        originParametersGuess.position(ctx.geoContext);
    m_vertexGuess = TVector3(originPositionGuess.x(), originPositionGuess.y(),
                             originPositionGuess.z());

    // Guessed bound track parameters
    Acts::BoundVector boundTrackParametersGuess =
        originParametersGuess.parameters();

    TArrayD boundTrackParsGuessData(Acts::eBoundSize);
    for (std::size_t i = 0; i < Acts::eBoundSize; i++) {
      boundTrackParsGuessData[i] = boundTrackParametersGuess(i);
    }
    m_boundTrackParametersGuess.Use(0, Acts::eBoundSize,
                                    boundTrackParsGuessData.GetArray());

    // Guessed bound errors
    Acts::BoundMatrix boundTrackCovGuess =
        originParametersGuess.covariance().value();
    TArrayD boundTrackCovGuessData(Acts::eBoundSize * Acts::eBoundSize);
    for (std::size_t i = 0; i < Acts::eBoundSize * Acts::eBoundSize; i++) {
      boundTrackCovGuessData[i] = boundTrackCovGuess(i);
    }
    m_boundTrackCovGuess.Use(Acts::eBoundSize, Acts::eBoundSize,
                             boundTrackCovGuessData.GetArray());

    // ----------------------------------------------
    // Estimated track parameters

    // Estimated origin momentum
    Acts::Vector3 originMomentumEst = track.momentum();
    m_originMomentumEst.SetPxPyPzE(
        originMomentumEst.x(), originMomentumEst.y(), originMomentumEst.z(),
        std::hypot(originMomentumEst.norm(), particleMass));

    // KF predicted vertex position
    Acts::Vector3 vertexEst = m_cfg.referenceSurface->localToGlobal(
        ctx.geoContext, {track.loc0(), track.loc1()},
        originMomentumEst.normalized());
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
    m_pdgId = originParametersGuess.particleHypothesis().absolutePdg();

    // Get charge
    m_charge = originParametersGuess.charge();

    // Get DoFs
    m_ndf = track.nDoF();

    // Get track ID
    m_trackId = trackIndices.trackId;

    // Iterate over the track states
    m_chi2Predicted = 0;
    m_chi2Filtered = 0;
    m_chi2Smoothed = 0;
    for (const auto& state : track.trackStatesReversed()) {
      // Skip the states without meaningful information
      if (!state.hasProjector()) {
        continue;
      }

      Acts::SourceLink stateSourceLink = state.getUncalibratedSourceLink();
      bool extended = (stateSourceLink.type() == typeid(ExtendedSourceLink));

      // Get the state reference surface
      const auto& referenceSurface = state.referenceSurface();

      // Get the source link index and flags
      std::size_t sourceLinkIdx =
          extended ? SimpleSourceLink(stateSourceLink.get<ExtendedSourceLink>())
                         .index()
                   : stateSourceLink.get<SimpleSourceLink>().index();
      bool backwards =
          extended ? stateSourceLink.get<ExtendedSourceLink>().isBackwards()
                   : false;

      m_geometryIds.push_back(referenceSurface.geometryId().sensitive());

      // ---------------------------------------------
      // Track hit info
      Acts::ActsDynamicVector calibratedPars = state.effectiveCalibrated();

      // Get the measurement hit
      Acts::Vector2 measHit = state.effectiveCalibrated().head(2);

      // Transform the hits to the global coordinates
      Acts::Vector3 measHitGlobal = referenceSurface.localToGlobal(
          ctx.geoContext, measHit, Acts::Vector3(1, 0, 0));

      // Get the measurement angles
      Acts::Vector2 measAngle = Acts::Vector2::Zero();
      if (extended) {
        measAngle = state.effectiveCalibrated().tail(2);
      }
      double measPhi = measAngle(0);
      double measTheta = measAngle(1);

      // Hit covariance
      Acts::ActsDynamicMatrix measCov = state.effectiveCalibratedCovariance();
      Acts::SquareMatrix2 measHitCov = measCov.topLeftCorner(2, 2);

      Acts::SquareMatrix2 measAngleCov = Acts::SquareMatrix2::Zero();
      if (extended) {
        measAngleCov = measCov.bottomRightCorner(2, 2);
      }

      // Store local measurement hits
      m_trackHitsLocal.emplace_back(measHit.x(), measHit.y());

      // Store global measurement hits
      m_trackHitsGlobal.emplace_back(measHitGlobal.x(), measHitGlobal.y(),
                                     measHitGlobal.z());

      // Store measurement direction
      m_onSurfaceTrackDirection.emplace_back(
          std::sin(measTheta) * std::cos(measPhi),
          std::sin(measTheta) * std::sin(measPhi), std::cos(measTheta));

      // Store measurement hit covariance
      TMatrixD trackHitCov(2, 2);
      TArrayD trackHitCovData(4);
      for (std::size_t i = 0; i < 4; i++) {
        trackHitCovData[i] = measHitCov(i);
      }
      trackHitCov.Use(2, 2, trackHitCovData.GetArray());
      m_trackHitCovs.push_back(trackHitCov);

      // Store measurement angle covariance
      TMatrixD trackAngleCov(2, 2);
      TArrayD trackAngleCovData(4);
      for (std::size_t i = 0; i < 4; i++) {
        trackAngleCovData[i] = measAngleCov(i);
      }
      trackAngleCov.Use(2, 2, trackAngleCovData.GetArray());
      m_trackAngleCovs.push_back(trackAngleCov);

      // ---------------------------------------------
      // Predicted state info
      Acts::BoundVector predictedPars = state.predicted();
      Acts::ActsDynamicMatrix effectiveProjector = state.effectiveProjector();

      // Get the predicted measurement hit
      Acts::ActsDynamicVector predictedMeas =
          effectiveProjector * predictedPars;
      Acts::Vector2 predictedHit = predictedMeas.head(2);

      // Transform the predicted hits to the global coordinates
      Acts::Vector3 predictedHitGlobal = referenceSurface.localToGlobal(
          ctx.geoContext, predictedHit, Acts::Vector3(1, 0, 0));

      // Get the predicted angles and momentum
      Acts::Vector2 predictedAngle = Acts::Vector2::Zero();
      if (extended) {
        predictedAngle = predictedMeas.tail(2);
      }
      double predictedOnSurfaceAbsMom =
          std::abs(1.0 / predictedPars(Acts::eBoundQOverP));

      // Get the residuals between the measurements and the predicted hits
      Acts::Vector2 predictedHitResidual = measHit - predictedHit;

      // Get the residuals between the measurements and the predicted angles
      Acts::Vector2 predictedAngleResidual = measAngle - predictedAngle;

      // Covariance with respect to the true track hits
      Acts::ActsDynamicMatrix predictedMeasCovTruth =
          effectiveProjector * state.predictedCovariance() *
          effectiveProjector.transpose();

      Acts::SquareMatrix2 predictedHitCovTruth =
          predictedMeasCovTruth.topLeftCorner(2, 2);
      Acts::SquareMatrix2 predictedAngleCovTruth = Acts::SquareMatrix2::Zero();
      if (extended) {
        predictedAngleCovTruth = predictedMeasCovTruth.bottomRightCorner(2, 2);
      }

      // Covariance with respect to the measurement hits
      Acts::ActsDynamicMatrix predictedMeasCov =
          measCov - predictedMeasCovTruth;

      Acts::SquareMatrix2 predictedHitCov =
          predictedMeasCov.topLeftCorner(2, 2);
      Acts::SquareMatrix2 predictedAngleCov = Acts::SquareMatrix2::Zero();
      if (extended) {
        predictedAngleCov = predictedMeasCov.bottomRightCorner(2, 2);
      }

      // Extract diagonals
      Acts::Vector2 predictedHitDiagTruth =
          predictedHitCovTruth.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
      Acts::Vector2 predictedAngleDiagTruth = Acts::Vector2::Zero();
      if (extended) {
        predictedAngleDiagTruth = predictedAngleCovTruth.cwiseAbs()
                                      .diagonal()
                                      .cwiseInverse()
                                      .cwiseSqrt();
      }

      Acts::Vector2 predictedHitDiag =
          predictedHitCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
      Acts::Vector2 predictedAngleDiag = Acts::Vector2::Zero();
      if (extended) {
        predictedAngleDiag =
            predictedAngleCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
      }

      // KF pulls with respect to the measurements
      Acts::Vector2 predictedHitPull =
          predictedHitDiag.cwiseProduct(measHit - predictedHit);

      Acts::Vector2 predictedAnglePull =
          predictedAngleDiag.cwiseProduct(measAngle - predictedAngle);

      // Store the KF predicted local hits
      m_predictedTrackHitsLocal.emplace_back(predictedHit.x(),
                                             predictedHit.y());

      // Store the KF predicted global hits
      m_predictedTrackHitsGlobal.emplace_back(predictedHitGlobal.x(),
                                              predictedHitGlobal.y(),
                                              predictedHitGlobal.z());

      // Store the KF predicted on surface momentum
      m_predictedOnSurfaceMomentum.emplace_back(
          predictedOnSurfaceAbsMom * std::sin(predictedAngle(1)) *
              std::cos(predictedAngle(0)),
          predictedOnSurfaceAbsMom * std::sin(predictedAngle(1)) *
              std::sin(predictedAngle(0)),
          predictedOnSurfaceAbsMom * std::cos(predictedAngle(1)),
          std::hypot(predictedOnSurfaceAbsMom, particleMass));

      // Store the residuals with respect to the measurements
      m_predictedHitResiduals.emplace_back(predictedHitResidual.x(),
                                           predictedHitResidual.y());

      m_predictedAngleResiduals.emplace_back(predictedAngleResidual(0),
                                             predictedAngleResidual(1));

      // Store the pulls with respect to the measurements
      m_predictedHitPulls.emplace_back(predictedHitPull.x(),
                                       predictedHitPull.y());

      m_predictedAnglePulls.emplace_back(predictedAnglePull(0),
                                         predictedAnglePull(1));

      // Add to the track chi2
      m_chi2Predicted += predictedHitPull.dot(predictedHitPull) +
                         predictedAnglePull.dot(predictedAnglePull);

      // ---------------------------------------------
      // Filtered state info
      if (state.hasFiltered()) {
        Acts::BoundVector filteredPars = state.filtered();
        Acts::ActsDynamicMatrix effectiveProjector = state.effectiveProjector();

        // Get the filtered measurement hit
        Acts::ActsDynamicVector filteredMeas =
            effectiveProjector * filteredPars;
        Acts::Vector2 filteredHit = filteredMeas.head(2);

        // Transform the filtered hits to the global coordinates
        Acts::Vector3 filteredHitGlobal = referenceSurface.localToGlobal(
            ctx.geoContext, filteredHit, Acts::Vector3(1, 0, 0));

        // Get the filtered angles and momentum
        Acts::Vector2 filteredAngle = Acts::Vector2::Zero();
        if (extended) {
          filteredAngle = filteredMeas.tail(2);
        }
        double filteredOnSurfaceAbsMom =
            std::abs(1.0 / filteredPars[Acts::eBoundQOverP]);

        // Get the residuals between the measurements and the filtered hits
        Acts::Vector2 filteredHitResidual = measHit - filteredHit;

        // Get the residuals between the measurements and the filtered angles
        Acts::Vector2 filteredAngleResidual = measAngle - filteredAngle;

        // Covariance with respect to the true track hits
        Acts::ActsDynamicMatrix filteredMeasCovTruth =
            effectiveProjector * state.filteredCovariance() *
            effectiveProjector.transpose();

        Acts::SquareMatrix2 filteredHitCovTruth =
            filteredMeasCovTruth.topLeftCorner(2, 2);
        Acts::SquareMatrix2 filteredAngleCovTruth = Acts::SquareMatrix2::Zero();
        if (extended) {
          filteredAngleCovTruth = filteredMeasCovTruth.bottomRightCorner(2, 2);
        }

        // Covariance with respect to the measurement hits
        Acts::ActsDynamicMatrix filteredMeasCov =
            measCov - filteredMeasCovTruth;

        Acts::SquareMatrix2 filteredHitCov =
            filteredMeasCov.topLeftCorner(2, 2);
        Acts::SquareMatrix2 filteredAngleCov = Acts::SquareMatrix2::Zero();
        if (extended) {
          filteredAngleCov = filteredMeasCov.bottomRightCorner(2, 2);
        }

        // Extract diagonals
        Acts::Vector2 filteredHitDiagTruth = filteredHitCovTruth.cwiseAbs()
                                                 .diagonal()
                                                 .cwiseInverse()
                                                 .cwiseSqrt();
        Acts::Vector2 filteredAngleDiagTruth = Acts::Vector2::Zero();
        if (extended) {
          filteredAngleDiagTruth = filteredAngleCovTruth.cwiseAbs()
                                       .diagonal()
                                       .cwiseInverse()
                                       .cwiseSqrt();
        }

        Acts::Vector2 filteredHitDiag =
            filteredHitCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
        Acts::Vector2 filteredAngleDiag = Acts::Vector2::Zero();
        if (extended) {
          filteredAngleDiag =
              filteredAngleCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
        }

        // KF pulls with respect to the measurements
        Acts::Vector2 filteredHitPull =
            filteredHitDiag.cwiseProduct(measHit - filteredHit);

        Acts::Vector2 filteredAnglePull =
            filteredAngleDiag.cwiseProduct(measAngle - filteredAngle);

        // Store the KF filtered local hits
        m_filteredTrackHitsLocal.emplace_back(filteredHit.x(), filteredHit.y());

        // Store the KF filtered global hits
        m_filteredTrackHitsGlobal.emplace_back(filteredHitGlobal.x(),
                                               filteredHitGlobal.y(),
                                               filteredHitGlobal.z());

        // Store the KF filtered on surface momentum
        m_filteredOnSurfaceMomentum.emplace_back(
            filteredOnSurfaceAbsMom * std::sin(filteredAngle(1)) *
                std::cos(filteredAngle(0)),
            filteredOnSurfaceAbsMom * std::sin(filteredAngle(1)) *
                std::sin(filteredAngle(0)),
            filteredOnSurfaceAbsMom * std::cos(filteredAngle(1)),
            std::hypot(filteredOnSurfaceAbsMom, particleMass));

        // Store the residuals with respect to the measurements
        m_filteredHitResiduals.emplace_back(filteredHitResidual.x(),
                                            filteredHitResidual.y());

        m_filteredAngleResiduals.emplace_back(filteredAngleResidual(0),
                                              filteredAngleResidual(1));

        // Store the pulls with respect to the measurements
        m_filteredHitPulls.emplace_back(filteredHitPull.x(),
                                        filteredHitPull.y());

        m_filteredAnglePulls.emplace_back(filteredAnglePull(0),
                                          filteredAnglePull(1));

        // Add to the track chi2
        m_chi2Filtered += filteredHitPull.dot(filteredHitPull) +
                          filteredAnglePull.dot(filteredAnglePull);
      }

      // ---------------------------------------------
      // Smoothed state info
      if (state.hasSmoothed()) {
        Acts::BoundVector smoothedPars = state.smoothed();
        Acts::ActsDynamicMatrix effectiveProjector = state.effectiveProjector();

        // Get the smoothed measurement hit
        Acts::ActsDynamicVector smoothedMeas =
            effectiveProjector * smoothedPars;
        Acts::Vector2 smoothedHit = smoothedMeas.head(2);

        // Transform the smoothed hits to the global coordinates
        Acts::Vector3 smoothedHitGlobal = referenceSurface.localToGlobal(
            ctx.geoContext, smoothedHit, Acts::Vector3(1, 0, 0));

        // Get the smoothed angles and momentum
        Acts::Vector2 smoothedAngle = Acts::Vector2::Zero();
        if (extended) {
          smoothedAngle = smoothedMeas.tail(2);
        }
        double smoothedOnSurfaceAbsMom =
            std::abs(1.0 / smoothedPars[Acts::eBoundQOverP]);

        // Get the residuals between the measurements and the smoothed hits
        Acts::Vector2 smoothedHitResidual = measHit - smoothedHit;

        // Get the residuals between the measurements and the smoothed angles
        Acts::Vector2 smoothedAngleResidual = measAngle - smoothedAngle;

        // Covariance with respect to the true track hits
        Acts::ActsDynamicMatrix smoothedMeasCovTruth =
            effectiveProjector * state.smoothedCovariance() *
            effectiveProjector.transpose();

        Acts::SquareMatrix2 smoothedHitCovTruth =
            smoothedMeasCovTruth.topLeftCorner(2, 2);
        Acts::SquareMatrix2 smoothedAngleCovTruth = Acts::SquareMatrix2::Zero();
        if (extended) {
          smoothedAngleCovTruth = smoothedMeasCovTruth.bottomRightCorner(2, 2);
        }

        // Covariance with respect to the measurement hits
        Acts::ActsDynamicMatrix smoothedMeasCov =
            measCov - smoothedMeasCovTruth;

        Acts::SquareMatrix2 smoothedHitCov =
            smoothedMeasCov.topLeftCorner(2, 2);
        Acts::SquareMatrix2 smoothedAngleCov = Acts::SquareMatrix2::Zero();
        if (extended) {
          smoothedAngleCov = smoothedMeasCov.bottomRightCorner(2, 2);
        }

        // Extract diagonals
        Acts::Vector2 smoothedHitDiagTruth = smoothedHitCovTruth.cwiseAbs()
                                                 .diagonal()
                                                 .cwiseInverse()
                                                 .cwiseSqrt();
        Acts::Vector2 smoothedAngleDiagTruth = Acts::Vector2::Zero();
        if (extended) {
          smoothedAngleDiagTruth = smoothedAngleCovTruth.cwiseAbs()
                                       .diagonal()
                                       .cwiseInverse()
                                       .cwiseSqrt();
        }

        Acts::Vector2 smoothedHitDiag =
            smoothedHitCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
        Acts::Vector2 smoothedAngleDiag = Acts::Vector2::Zero();
        if (extended) {
          smoothedAngleDiag =
              smoothedAngleCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
        }

        // KF pulls with respect to the measurements
        Acts::Vector2 smoothedHitPull =
            smoothedHitDiag.cwiseProduct(measHit - smoothedHit);

        Acts::Vector2 smoothedAnglePull =
            smoothedAngleDiag.cwiseProduct(measAngle - smoothedAngle);

        // Store the KF smoothed local hits
        m_smoothedTrackHitsLocal.emplace_back(smoothedHit.x(), smoothedHit.y());

        // Store the KF smoothed global hits
        m_smoothedTrackHitsGlobal.emplace_back(smoothedHitGlobal.x(),
                                               smoothedHitGlobal.y(),
                                               smoothedHitGlobal.z());

        // Store the KF smoothed on surface momentum
        m_smoothedOnSurfaceMomentum.emplace_back(
            smoothedOnSurfaceAbsMom * std::sin(smoothedAngle(1)) *
                std::cos(smoothedAngle(0)),
            smoothedOnSurfaceAbsMom * std::sin(smoothedAngle(1)) *
                std::sin(smoothedAngle(0)),
            smoothedOnSurfaceAbsMom * std::cos(smoothedAngle(1)),
            std::hypot(smoothedOnSurfaceAbsMom, particleMass));

        // Store the residuals with respect to the measurements
        m_smoothedHitResiduals.emplace_back(smoothedHitResidual.x(),
                                            smoothedHitResidual.y());

        m_smoothedAngleResiduals.emplace_back(smoothedAngleResidual(0),
                                              smoothedAngleResidual(1));

        // Store the pulls with respect to the measurements
        m_smoothedHitPulls.emplace_back(smoothedHitPull.x(),
                                        smoothedHitPull.y());

        m_smoothedAnglePulls.emplace_back(smoothedAnglePull(0),
                                          smoothedAnglePull(1));

        // Add to the track chi2
        m_chi2Smoothed += smoothedHitPull.dot(smoothedHitPull) +
                          smoothedAnglePull.dot(smoothedAnglePull);
      }
    }

    // Fill the tree
    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
