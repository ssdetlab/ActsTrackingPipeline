#include "TrackingPipeline/Io/E320RootSimTrackWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <stdexcept>

#include "TrackingPipeline/EventData/ExtendedSourceLink.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/E320GeometryOptions.hpp"
#include "TrackingPipeline/MagneticField/ConstantMagField.hpp"
#include "TrackingPipeline/MagneticField/IdealQuadrupoleMagField.hpp"
#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"

E320::E320RootSimTrackWriter::E320RootSimTrackWriter(const Config& config,
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

  // Magnet configuration as seen by track
  m_tree->Branch("quad1Grad", &m_quad1Grad, bufSize, splitLvl);
  m_tree->Branch("quad2Grad", &m_quad2Grad, bufSize, splitLvl);
  m_tree->Branch("quad3Grad", &m_quad3Grad, bufSize, splitLvl);
  m_tree->Branch("xCorrectorStrength", &m_xCorrectorStrength, bufSize,
                 splitLvl);
  m_tree->Branch("dipoleStrength", &m_dipoleStrength, bufSize, splitLvl);

  // True track hits in the surface frame
  m_tree->Branch("trueTrackHitsLocal", &m_trueTrackHitsLocal, bufSize,
                 splitLvl);

  // True track hits in the global frame
  m_tree->Branch("trueTrackHitsGlobal", &m_trueTrackHitsGlobal, bufSize,
                 splitLvl);

  // True track momentum on the surface in the track frame
  m_tree->Branch("onSurfaceMomentumTruth", &m_onSurfaceMomentumTruth, bufSize,
                 splitLvl);

  // Signal/background flag of the track's clusters
  m_tree->Branch("isSignal", &m_isSignal, bufSize, splitLvl);

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

  // KF residuals with respect to the true hits
  m_tree->Branch("truePredictedHitResiduals", &m_truePredictedHitResiduals,
                 bufSize, splitLvl);
  m_tree->Branch("trueFilteredHitResiduals", &m_trueFilteredHitResiduals,
                 bufSize, splitLvl);
  m_tree->Branch("trueSmoothedHitResiduals", &m_trueSmoothedHitResiduals,
                 bufSize, splitLvl);

  m_tree->Branch("truePredictedAngleResiduals", &m_truePredictedAngleResiduals,
                 bufSize, splitLvl);
  m_tree->Branch("trueFilteredAngleResiduals", &m_trueFilteredAngleResiduals,
                 bufSize, splitLvl);
  m_tree->Branch("trueSmoothedAngleResiduals", &m_trueSmoothedAngleResiduals,
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

  // KF pulls with respect to the true hits
  m_tree->Branch("truePredictedHitPulls", &m_truePredictedHitPulls, bufSize,
                 splitLvl);
  m_tree->Branch("trueFilteredHitPulls", &m_trueFilteredHitPulls, bufSize,
                 splitLvl);
  m_tree->Branch("trueSmoothedHitPulls", &m_trueSmoothedHitPulls, bufSize,
                 splitLvl);

  m_tree->Branch("truePredictedAnglePulls", &m_truePredictedAnglePulls, bufSize,
                 splitLvl);
  m_tree->Branch("trueFilteredAnglePulls", &m_trueFilteredAnglePulls, bufSize,
                 splitLvl);
  m_tree->Branch("trueSmoothedAnglePulls", &m_trueSmoothedAnglePulls, bufSize,
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

  // True bound track parameters
  m_tree->Branch("boundTrackParametersTruth", &m_boundTrackParametersTruth,
                 bufSize, splitLvl);
  m_tree->Branch("boundTrackCovTruth", &m_boundTrackCovTruth, bufSize,
                 splitLvl);

  // Initial guess of the momentum at the IP
  m_tree->Branch("originMomentumGuess", &m_originMomentumGuess, bufSize,
                 splitLvl);

  // Initial guess of the vertex at the IP
  m_tree->Branch("vertexGuess", &m_vertexGuess, bufSize, splitLvl);

  // KF predicted momentum at the IP
  m_tree->Branch("originMomentumEst", &m_originMomentumEst, bufSize, splitLvl);

  // KF predicted vertex at the IP
  m_tree->Branch("vertexEst", &m_vertexEst, bufSize, splitLvl);

  // True momentum at the IP
  m_tree->Branch("originMomentumTruth", &m_originMomentumTruth, bufSize,
                 splitLvl);

  // True vertex at the IP
  m_tree->Branch("vertexTruth", &m_vertexTruth, bufSize, splitLvl);

  // Chi2 of the track with respect ot the measurement
  m_tree->Branch("chi2Predicted", &m_chi2Predicted, bufSize, splitLvl);
  m_tree->Branch("chi2Filtered", &m_chi2Filtered, bufSize, splitLvl);
  m_tree->Branch("chi2Smoothed", &m_chi2Smoothed, bufSize, splitLvl);

  // Number of degrees of freedom of the track
  m_tree->Branch("ndf", &m_ndf, bufSize, splitLvl);

  // Track IDs of track states' sim clusters
  m_tree->Branch("stateTrackId", &m_stateTrackId, bufSize, splitLvl);
  m_tree->Branch("stateParentTrackId", &m_stateParentTrackId, bufSize,
                 splitLvl);
  m_tree->Branch("stateRunId", &m_stateRunId, bufSize, splitLvl);

  // Track IDs of the track selected as the majority among the clusters' track
  // IDs
  m_tree->Branch("trackId", &m_trackId, bufSize, splitLvl);
  m_tree->Branch("parentTrackId", &m_parentTrackId, bufSize, splitLvl);
  m_tree->Branch("runId", &m_runId, bufSize, splitLvl);

  // Event ID
  m_tree->Branch("eventId", &m_eventId, bufSize, splitLvl);

  // True number of clusters in the track
  m_tree->Branch("trueTrackSize", &m_trueTrackSize, bufSize, splitLvl);

  // Number of captured clusters in the track
  m_tree->Branch("capturedTrackSize", &m_capturedTrackSize, bufSize, splitLvl);

  // PDG ID
  m_tree->Branch("pdgId", &m_pdgId, bufSize, splitLvl);

  // Charge
  m_tree->Branch("charge", &m_charge, bufSize, splitLvl);

  //------------------------------------------------------------------
  // Initialize the data handles
  m_inputTrackContainer.initialize(m_cfg.inputTrackContainer);
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputTrackParametersGuesses.initialize(m_cfg.inputTrackParametersGuesses);
  m_inputSimClusters.initialize(m_cfg.inputSimClusters);
}

ProcessCode E320::E320RootSimTrackWriter::finalize() {
  if (m_file != nullptr) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode E320::E320RootSimTrackWriter::write(const AlgorithmContext& ctx) {
  const auto& inputTrackContainer = m_inputTrackContainer(ctx);
  const auto& inputTracks = m_inputTracks(ctx);
  const auto& inputTrackParametersGuesses = m_inputTrackParametersGuesses(ctx);
  const auto& inputTruthClusters = m_inputSimClusters(ctx);

  std::lock_guard<std::mutex> lock(m_mutex);

  // Collect true track statistics
  auto trueTrackIds =
      inputTruthClusters |
      std::views::filter([](const auto& cl) { return cl.isSignal; }) |
      std::views::transform([](const auto& cl) { return cl.truthHits; }) |
      std::views::join | std::views::transform([](const auto& hit) -> TrackID {
        return {hit.trackId, hit.parentTrackId, hit.runId};
      });

  // Event meta data
  const auto& mctx = ctx.magFieldContext;
  const auto& goInst = *E320::GeometryOptions::instance();

  m_quad1Grad =
      goInst.quad1Gradient * Acts::UnitConstants::m / Acts::UnitConstants::T;
  m_quad2Grad =
      goInst.quad2Gradient * Acts::UnitConstants::m / Acts::UnitConstants::T;
  m_quad3Grad =
      goInst.quad3Gradient * Acts::UnitConstants::m / Acts::UnitConstants::T;

  // XCOR strength T
  m_xCorrectorStrength =
      goInst.xCorrectorFieldStrength / Acts::UnitConstants::T;

  // Dipole strength T
  m_dipoleStrength = goInst.dipoleFieldStrength / Acts::UnitConstants::T;

  if (mctx.hasValue()) {
    const auto& store = mctx.get<std::shared_ptr<MagneticFieldStore>&>();

    if (store->store.contains(goInst.quad1Id)) {
      m_quad1Grad = store->store.at(goInst.quad1Id)
                        .as<IdealQuadrupoleMagField::Cache>()
                        .m_gradient *
                    Acts::UnitConstants::m / Acts::UnitConstants::T;
    }
    if (store->store.contains(goInst.quad2Id)) {
      m_quad2Grad = store->store.at(goInst.quad2Id)
                        .as<IdealQuadrupoleMagField::Cache>()
                        .m_gradient *
                    Acts::UnitConstants::m / Acts::UnitConstants::T;
    }
    if (store->store.contains(goInst.quad3Id)) {
      m_quad3Grad = store->store.at(goInst.quad3Id)
                        .as<IdealQuadrupoleMagField::Cache>()
                        .m_gradient *
                    Acts::UnitConstants::m / Acts::UnitConstants::T;
    }
    if (store->store.contains(goInst.xCorrectorId)) {
      m_xCorrectorStrength = store->store.at(goInst.xCorrectorId)
                                 .as<ConstantMagField::Cache>()
                                 .m_field(goInst.xCorrectorDirIdx) /
                             Acts::UnitConstants::T;
    }
    if (store->store.contains(goInst.dipoleId)) {
      m_dipoleStrength = store->store.at(goInst.dipoleId)
                             .as<ConstantMagField::Cache>()
                             .m_field(goInst.dipoleDirIdx) /
                         Acts::UnitConstants::T;
    }
  }

  // Event ID
  m_eventId = ctx.eventNumber;

  // Iterate over the fitted tracks
  for (std::size_t idx = 0; idx < inputTracks.size(); idx++) {
    // Get track indices
    const auto& trackIndices = inputTracks.at(idx);

    // Get the track object and the track id
    const auto& track = inputTrackContainer.getTrack(trackIndices.trackIndex);

    std::size_t nStates = track.nTrackStates();

    // True track hits in the surface frame
    m_trueTrackHitsLocal.clear();
    m_trueTrackHitsLocal.reserve(nStates);

    // True track hits in the global frame
    m_trueTrackHitsGlobal.clear();
    m_trueTrackHitsGlobal.reserve(nStates);

    // True track momentum on the surface in the track frame
    m_onSurfaceMomentumTruth.clear();
    m_onSurfaceMomentumTruth.reserve(nStates);

    // Signal/background flag of the track's clusters
    m_isSignal.clear();
    m_isSignal.reserve(nStates);

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

    // KF residuals with respect to the true hits
    m_truePredictedHitResiduals.clear();
    m_truePredictedHitResiduals.reserve(nStates);

    m_trueFilteredHitResiduals.clear();
    m_trueFilteredHitResiduals.reserve(nStates);

    m_trueSmoothedHitResiduals.clear();
    m_trueSmoothedHitResiduals.reserve(nStates);

    m_truePredictedAngleResiduals.clear();
    m_truePredictedAngleResiduals.reserve(nStates);

    m_trueFilteredAngleResiduals.clear();
    m_trueFilteredAngleResiduals.reserve(nStates);

    m_trueSmoothedAngleResiduals.clear();
    m_trueSmoothedAngleResiduals.reserve(nStates);

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

    // KF pulls with respect to the true hits
    m_truePredictedHitPulls.clear();
    m_truePredictedHitPulls.reserve(nStates);

    m_trueFilteredHitPulls.clear();
    m_trueFilteredHitPulls.reserve(nStates);

    m_trueSmoothedHitPulls.clear();
    m_trueSmoothedHitPulls.reserve(nStates);

    m_truePredictedAnglePulls.clear();
    m_truePredictedAnglePulls.reserve(nStates);

    m_trueFilteredAnglePulls.clear();
    m_trueFilteredAnglePulls.reserve(nStates);

    m_trueSmoothedAnglePulls.clear();
    m_trueSmoothedAnglePulls.reserve(nStates);

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

    // Track IDs of track states' sim clusters
    m_stateTrackId.clear();
    m_stateTrackId.reserve(nStates);

    m_stateParentTrackId.clear();
    m_stateParentTrackId.reserve(nStates);

    m_stateRunId.clear();
    m_stateRunId.reserve(nStates);

    // ----------------------------------------------
    // Guess track parameters

    const auto& originParametersGuess =
        inputTrackParametersGuesses.at(trackIndices.originParametersGuessIndex);

    // Guessed IP momentum
    const auto& originMomentumGuess = originParametersGuess.momentum();
    double particleMass = originParametersGuess.particleHypothesis().mass();
    m_originMomentumGuess.SetPxPyPzE(
        originMomentumGuess.x(), originMomentumGuess.y(),
        originMomentumGuess.z(),
        std::hypot(originMomentumGuess.norm(), particleMass));

    // Guessed vertex
    const auto& ipPositionGuess =
        originParametersGuess.position(ctx.geoContext);
    m_vertexGuess =
        TVector3(ipPositionGuess.x(), ipPositionGuess.y(), ipPositionGuess.z());

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

    // Estimated IP momentum
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

    // Iterate over the track states
    std::map<TrackID, std::vector<int>> trackStateIds;
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

      const auto& cluster = inputTruthClusters.at(sourceLinkIdx);
      m_isSignal.push_back(static_cast<int>(cluster.isSignal));

      // Get the true hit info
      Acts::ActsDynamicVector calibratedPars = state.effectiveCalibrated();

      TrackID currentTrackId;
      Acts::BoundVector onSurfaceTruthParameters;
      if (cluster.truthHits.size() == 0 || !cluster.isSignal) {
        // No signal tracks in the cluster
        currentTrackId = std::make_tuple(-1, -1, -1);

        onSurfaceTruthParameters = Acts::BoundVector::Zero();
        onSurfaceTruthParameters.head(2) = calibratedPars.head(2);
        onSurfaceTruthParameters(Acts::eBoundQOverP) = -1;
      } else {
        // Found signal track in the cluster
        auto sig = std::ranges::find_if(cluster.truthHits, [](const auto& hit) {
          return (hit.trackId == 1);
        });
        onSurfaceTruthParameters = sig->truthParameters;
        currentTrackId =
            std::make_tuple(sig->trackId, sig->parentTrackId, sig->runId);
      }
      trackStateIds[currentTrackId].push_back(sourceLinkIdx);

      // Get true local track hit
      Acts::Vector2 trueTrackHit = onSurfaceTruthParameters.head(2);

      // Transform the hits to the global coordinates
      Acts::Vector3 trueHitGlobal =
          m_cfg.surfaceAccessor(cluster.sourceLink)
              ->localToGlobal(ctx.geoContext, trueTrackHit,
                              Acts::Vector3(1, 0, 0));

      // Get true local track angles
      Acts::Vector2 trueTrackAngle = onSurfaceTruthParameters.segment(2, 2);

      // ---------------------------------------------
      // Track hit info

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

      // Get the residuals between the true and the predicted hits
      Acts::Vector2 truePredictedHitResidual = trueTrackHit - predictedHit;

      // Get the residuals between the true and the predicted angles
      if (m_cfg.applyPhiCorrection) {
        double measPhi = trueTrackAngle(0);
        double trackPhi = predictedPars(Acts::eBoundPhi);
        double phiDiff = measPhi - trackPhi;
        if (phiDiff < 0 && std::abs(phiDiff) > M_PI) {
          trueTrackAngle(0) += 2 * M_PI;
        }
        if (phiDiff > M_PI) {
          trueTrackAngle(0) -= 2 * M_PI;
        }
      }
      Acts::Vector2 truePredictedAngleResidual =
          trueTrackAngle - predictedAngle;

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

      // KF pulls with respect to the true hits
      Acts::Vector2 truePredictedHitPull =
          predictedHitDiagTruth.cwiseProduct(trueTrackHit - predictedHit);

      Acts::Vector2 truePredictedAnglePull =
          predictedAngleDiagTruth.cwiseProduct(trueTrackAngle - predictedAngle);

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

      // Store the residuals with respect to the true hits
      m_truePredictedHitResiduals.emplace_back(truePredictedHitResidual.x(),
                                               truePredictedHitResidual.y());

      m_truePredictedAngleResiduals.emplace_back(truePredictedAngleResidual(0),
                                                 truePredictedAngleResidual(1));

      // Store the residuals with respect to the measurements
      m_predictedHitResiduals.emplace_back(predictedHitResidual.x(),
                                           predictedHitResidual.y());

      m_predictedAngleResiduals.emplace_back(predictedAngleResidual(0),
                                             predictedAngleResidual(1));

      // Store the pulls with respect to the true hits
      m_truePredictedHitPulls.emplace_back(truePredictedHitPull.x(),
                                           truePredictedHitPull.y());

      m_truePredictedAnglePulls.emplace_back(truePredictedAnglePull(0),
                                             truePredictedAnglePull(1));

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

        // Get the residuals between the true and the filtered hits
        Acts::Vector2 trueFilteredHitResidual = trueTrackHit - filteredHit;

        // Get the residuals between the true and the filtered angles
        Acts::Vector2 trueFilteredAngleResidual =
            trueTrackAngle - filteredAngle;

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

        // KF pulls with respect to the true hits
        Acts::Vector2 trueFilteredHitPull =
            filteredHitDiagTruth.cwiseProduct(trueTrackHit - filteredHit);

        Acts::Vector2 trueFilteredAnglePull =
            filteredAngleDiagTruth.cwiseProduct(trueTrackAngle - filteredAngle);

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

        // Store the residuals with respect to the true hits
        m_trueFilteredHitResiduals.emplace_back(trueFilteredHitResidual.x(),
                                                trueFilteredHitResidual.y());

        m_trueFilteredAngleResiduals.emplace_back(trueFilteredAngleResidual(0),
                                                  trueFilteredAngleResidual(1));

        // Store the residuals with respect to the measurements
        m_filteredHitResiduals.emplace_back(filteredHitResidual.x(),
                                            filteredHitResidual.y());

        m_filteredAngleResiduals.emplace_back(filteredAngleResidual(0),
                                              filteredAngleResidual(1));

        // Store the pulls with respect to the true hits
        m_trueFilteredHitPulls.emplace_back(trueFilteredHitPull.x(),
                                            trueFilteredHitPull.y());

        m_trueFilteredAnglePulls.emplace_back(trueFilteredAnglePull(0),
                                              trueFilteredAnglePull(1));

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

        // Get the residuals between the true and the smoothed hits
        Acts::Vector2 trueSmoothedHitResidual = trueTrackHit - smoothedHit;

        // Get the residuals between the true and the smoothed angles
        Acts::Vector2 trueSmoothedAngleResidual =
            trueTrackAngle - smoothedAngle;

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

        // KF pulls with respect to the true hits
        Acts::Vector2 trueSmoothedHitPull =
            smoothedHitDiagTruth.cwiseProduct(trueTrackHit - smoothedHit);

        Acts::Vector2 trueSmoothedAnglePull =
            smoothedAngleDiagTruth.cwiseProduct(trueTrackAngle - smoothedAngle);

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

        // Store the residuals with respect to the true hits
        m_trueSmoothedHitResiduals.emplace_back(trueSmoothedHitResidual.x(),
                                                trueSmoothedHitResidual.y());

        m_trueSmoothedAngleResiduals.emplace_back(trueSmoothedAngleResidual(0),
                                                  trueSmoothedAngleResidual(1));

        // Store the residuals with respect to the measurements
        m_smoothedHitResiduals.emplace_back(smoothedHitResidual.x(),
                                            smoothedHitResidual.y());

        m_smoothedAngleResiduals.emplace_back(smoothedAngleResidual(0),
                                              smoothedAngleResidual(1));

        // Store the pulls with respect to the true hits
        m_trueSmoothedHitPulls.emplace_back(trueSmoothedHitPull.x(),
                                            trueSmoothedHitPull.y());

        m_trueSmoothedAnglePulls.emplace_back(trueSmoothedAnglePull(0),
                                              trueSmoothedAnglePull(1));

        // Store the pulls with respect to the measurements
        m_smoothedHitPulls.emplace_back(smoothedHitPull.x(),
                                        smoothedHitPull.y());

        m_smoothedAnglePulls.emplace_back(smoothedAnglePull(0),
                                          smoothedAnglePull(1));

        // Add to the track chi2
        m_chi2Smoothed += smoothedHitPull.dot(smoothedHitPull) +
                          smoothedAnglePull.dot(smoothedAnglePull);
      }

      // Store the true local hits
      m_trueTrackHitsLocal.emplace_back(trueTrackHit.x(), trueTrackHit.y());

      // Store the true global hits
      m_trueTrackHitsGlobal.emplace_back(trueHitGlobal.x(), trueHitGlobal.y(),
                                         trueHitGlobal.z());

      // Store true on surface momentum
      TLorentzVector onSurfaceMomTruth;
      double trueOnSurfaceAbsMom =
          std::abs(1. / onSurfaceTruthParameters(Acts::eBoundQOverP));
      onSurfaceMomTruth.SetPxPyPzE(
          trueOnSurfaceAbsMom * std::sin(trueTrackAngle(1)) *
              std::cos(trueTrackAngle(0)),
          trueOnSurfaceAbsMom * std::sin(trueTrackAngle(1)) *
              std::sin(trueTrackAngle(0)),
          trueOnSurfaceAbsMom * std::cos(trueTrackAngle(1)),
          std::hypot(trueOnSurfaceAbsMom, particleMass));
      m_onSurfaceMomentumTruth.emplace_back(onSurfaceMomTruth);

      // Store state track ID
      m_stateTrackId.push_back(std::get<0>(currentTrackId));
      m_stateParentTrackId.push_back(std::get<1>(currentTrackId));
      m_stateRunId.push_back(std::get<2>(currentTrackId));
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
      TArrayD boundTrackParsTruthData(Acts::eBoundSize);
      for (std::size_t i = 0; i < Acts::eBoundSize; i++) {
        boundTrackParsTruthData[i] = 0;
      }
      TArrayD boundTrackCovTruthData(Acts::eBoundSize * Acts::eBoundSize);
      for (std::size_t i = 0; i < Acts::eBoundSize * Acts::eBoundSize; i++) {
        boundTrackCovTruthData[i] = 0;
      }

      m_originMomentumTruth.SetPxPyPzE(0, 0, 0, 0);
      m_vertexTruth.SetXYZ(0, 0, 0);

      m_boundTrackParametersTruth.Use(0, Acts::eBoundSize,
                                      boundTrackParsTruthData.GetArray());
      m_boundTrackCovTruth.Use(Acts::eBoundSize, Acts::eBoundSize,
                               boundTrackCovTruthData.GetArray());

      m_trueTrackSize = 0;
      m_capturedTrackSize = 0;

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

      const auto& ipParametersTruth = pivotHit->ipParameters;

      // Truth IP momentum
      const auto& originMomentumTruth = ipParametersTruth.momentum();
      double particleMass = ipParametersTruth.particleHypothesis().mass();
      m_originMomentumTruth.SetPxPyPzE(
          originMomentumTruth.x(), originMomentumTruth.y(),
          originMomentumTruth.z(),
          std::hypot(originMomentumTruth.norm(), particleMass));

      // Truth vertex
      const auto& ipPositionTruth = ipParametersTruth.position(ctx.geoContext);
      m_vertexTruth = TVector3(ipPositionTruth.x(), ipPositionTruth.y(),
                               ipPositionTruth.z());

      // Truth bound track parameters
      Acts::BoundVector boundTrackParametersTruth =
          ipParametersTruth.parameters();

      TArrayD boundTrackParsTruthData(Acts::eBoundSize);
      for (std::size_t i = 0; i < Acts::eBoundSize; i++) {
        boundTrackParsTruthData[i] = boundTrackParametersTruth(i);
      }
      m_boundTrackParametersTruth.Use(0, Acts::eBoundSize,
                                      boundTrackParsTruthData.GetArray());

      // Truth bound errors
      Acts::BoundMatrix boundTrackCovTruth =
          ipParametersTruth.covariance().value();
      TArrayD boundTrackCovTruthData(Acts::eBoundSize * Acts::eBoundSize);
      for (std::size_t i = 0; i < Acts::eBoundSize * Acts::eBoundSize; i++) {
        boundTrackCovTruthData[i] = boundTrackCovTruth(i);
      }
      m_boundTrackCovTruth.Use(Acts::eBoundSize, Acts::eBoundSize,
                               boundTrackCovTruthData.GetArray());

      // Compute matching degree
      m_trueTrackSize = std::ranges::count(trueTrackIds, refTrackId->first);
      m_capturedTrackSize = refTrackId->second.size();

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
