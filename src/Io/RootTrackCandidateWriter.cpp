#include "TrackingPipeline/Io/RootTrackCandidateWriter.hpp"

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

RootTrackCandidateWriter::RootTrackCandidateWriter(const Config& config,
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
  m_tree->Branch("trackHits", &m_trackHits, buf_size, split_lvl);

  // CKF predicted track hits
  m_tree->Branch("predictedTrackHits", &m_predictedTrackHits, buf_size,
                 split_lvl);
  m_tree->Branch("filteredTrackHits", &m_filteredTrackHits, buf_size,
                 split_lvl);

  // CKF residuals with respect to the measurements
  m_tree->Branch("predictedResiduals", &m_predictedResiduals, buf_size,
                 split_lvl);
  m_tree->Branch("filteredResiduals", &m_filteredResiduals, buf_size,
                 split_lvl);

  // CKF pulls with respect to the measurements
  m_tree->Branch("predictedPulls", &m_predictedPulls, buf_size, split_lvl);
  m_tree->Branch("filteredPulls", &m_filteredPulls, buf_size, split_lvl);

  // Chi2 and ndf of the fitted track
  m_tree->Branch("chi2Predicted", &m_chi2Predicted, buf_size, split_lvl);
  m_tree->Branch("chi2Filtered", &m_chi2Filtered, buf_size, split_lvl);
  m_tree->Branch("ndf", &m_ndf, "ndf/I");

  // Track ID
  m_tree->Branch("trackId", &m_trackId, "trackId/I");

  // Event ID
  m_tree->Branch("eventId", &m_eventId, "eventId/I");

  //------------------------------------------------------------------
  // Initialize the data handles
  m_inputTrackCandidates.initialize(m_cfg.inputTrackCandidates);
}

ProcessCode RootTrackCandidateWriter::finalize() {
  if (m_file) {
    m_file->Write();
    m_file->Close();
    delete m_file;
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootTrackCandidateWriter::write(const AlgorithmContext& ctx) {
  auto inputCandidates = m_inputTrackCandidates(ctx);

  std::lock_guard<std::mutex> lock(m_mutex);

  m_eventId = ctx.eventNumber;

  // Iterate over the fitted tracks
  for (std::size_t tid = 0; tid < inputCandidates.size(); tid++) {
    // Get the track object and the track id
    const auto& track = inputCandidates.getTrack(tid);

    // Track hits from the measurements
    std::vector<TVector3> trackHits;

    // CKF predicted track hits
    std::vector<TVector3> predictedTrackHits;
    std::vector<TVector3> filteredTrackHits;

    // CKF residuals with respect to the measurements
    std::vector<TVector3> predictedResiduals;
    std::vector<TVector3> filteredResiduals;

    // CKF pulls with respect to the measurements
    std::vector<TVector3> predictedPulls;
    std::vector<TVector3> filteredPulls;

    // Chi2 of different track states
    std::vector<double> chi2Predicted;
    std::vector<double> chi2Filtered;

    // Iterate over the track states
    for (const auto& state : track.trackStatesReversed()) {
      // Skip the states without meaningful information
      if (!state.hasProjector()) {
        continue;
      }

      // Get the measurements source link
      auto sl = state.getUncalibratedSourceLink();
      auto ssl = sl.get<SimpleSourceLink>();

      // Get the measurements hit
      auto hit = state.effectiveCalibrated();

      // Project onto the prediction space
      auto predictedHit = state.effectiveProjector() * state.predicted();
      auto filteredHit = state.effectiveProjector() * state.filtered();

      // Transform the hits to the global coordinates
      auto hitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, hit, Acts::Vector3(1, 0, 0));
      auto predictedHitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, predictedHit, Acts::Vector3(1, 0, 0));
      auto filteredHitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, filteredHit, Acts::Vector3(1, 0, 0));

      // Get the residuals between the measurements and the predicted hits
      auto predictedResidual = hitGlobal - predictedHitGlobal;
      auto filteredResidual = hitGlobal - filteredHitGlobal;

      // CKF predicted covariances
      auto measurementCov = state.effectiveCalibratedCovariance();

      // With respect to measurement
      auto predictedCov = state.effectiveProjector() *
                          state.predictedCovariance() *
                          state.effectiveProjector().transpose();

      auto filteredCov = state.effectiveProjector() *
                             state.filteredCovariance() *
                             state.effectiveProjector().transpose() -
                         measurementCov;

      // Extract diagonals
      auto predictedDiag =
          predictedCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
      auto filteredDiag =
          filteredCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

      // CKF pulls with respect to the measurements
      auto predictedPull = predictedDiag.cwiseProduct(hit - predictedHit);
      auto filteredPull = filteredDiag.cwiseProduct(hit - filteredHit);

      // Store the measurements hits
      trackHits.push_back(
          TVector3(hitGlobal.x(), hitGlobal.y(), hitGlobal.z()));

      // Store the KF predicted hits
      predictedTrackHits.push_back(TVector3(predictedHitGlobal.x(),
                                            predictedHitGlobal.y(),
                                            predictedHitGlobal.z()));
      filteredTrackHits.push_back(TVector3(
          filteredHitGlobal.x(), filteredHitGlobal.y(), filteredHitGlobal.z()));

      // Store the residuals with respect to the measurements
      predictedResiduals.push_back(TVector3(
          predictedResidual.x(), predictedResidual.y(), predictedResidual.z()));
      filteredResiduals.push_back(TVector3(
          filteredResidual.x(), filteredResidual.y(), filteredResidual.z()));

      // Store the pulls with respect to the measurements
      predictedPulls.push_back(
          TVector3(predictedPull.x(), 0, -predictedPull.y()));
      filteredPulls.push_back(TVector3(filteredPull.x(), 0, -filteredPull.y()));

      // Chi2 of the track state
      chi2Predicted.push_back(predictedPull.dot(predictedPull));
      chi2Filtered.push_back(filteredPull.dot(filteredPull));
    }
    double me = 0.511 * Acts::UnitConstants::MeV;

    // Measurement hits
    m_trackHits = trackHits;

    // CKF predicted track hits
    m_predictedTrackHits = predictedTrackHits;
    m_filteredTrackHits = filteredTrackHits;

    // CKF residuals with respect to the measurements
    m_predictedResiduals = predictedResiduals;
    m_filteredResiduals = filteredResiduals;

    // CKF pulls with respect to the measurements
    m_predictedPulls = predictedPulls;
    m_filteredPulls = filteredPulls;

    // Chi2 of the track
    // with respect ot the
    // measurement
    m_chi2Predicted = chi2Predicted;
    m_chi2Filtered = chi2Filtered;

    // Number of degrees of freedom
    m_ndf = track.nDoF();

    // Track Id
    m_trackId = tid;

    // Fill the tree
    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
