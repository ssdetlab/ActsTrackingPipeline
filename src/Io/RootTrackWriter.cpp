#include "TrackingPipeline/Io/RootTrackWriter.hpp"

#include <algorithm>
#include <ranges>

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
  m_tree->Branch("trackHits", &m_trackHits, buf_size, split_lvl);

  // KF predicted track hits
  m_tree->Branch("predictedTrackHits", &m_predictedTrackHits, buf_size,
                 split_lvl);
  m_tree->Branch("filteredTrackHits", &m_filteredTrackHits, buf_size,
                 split_lvl);
  m_tree->Branch("smoothedTrackHits", &m_smoothedTrackHits, buf_size,
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

  // KF predicted momentum at the IP
  m_tree->Branch("ipMomentum", &m_ipMomentum);
  m_tree->Branch("ipMomentumError", &m_ipMomentumError);
  m_tree->Branch("vertex", &m_vertex);
  m_tree->Branch("vertexError", &m_vertexError);

  // Chi2 and ndf of the fitted track
  m_tree->Branch("chi2Predicted", &m_chi2Predicted, "chi2Predicted/D");
  m_tree->Branch("chi2Filtered", &m_chi2Filtered, "chi2Filtered/D");
  m_tree->Branch("chi2Smoothed", &m_chi2Smoothed, "chi2Smoothed/D");
  m_tree->Branch("ndf", &m_ndf, "ndf/I");

  // Track ID
  m_tree->Branch("trackId", &m_trackId, "trackId/I");

  // Event ID
  m_tree->Branch("eventId", &m_eventId, "eventId/I");

  //------------------------------------------------------------------
  // Initialize the data handles
  m_inputTracks.initialize(m_cfg.inputTracks);
}

ProcessCode RootTrackWriter::finalize() {
  if (m_file) {
    m_file->Write();
    m_file->Close();
    delete m_file;
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootTrackWriter::write(const AlgorithmContext& ctx) {
  auto inputTracks = m_inputTracks(ctx);

  std::lock_guard<std::mutex> lock(m_mutex);

  m_eventId = ctx.eventNumber;

  // Iterate over the fitted tracks
  for (std::size_t tid = 0; tid < inputTracks.size(); tid++) {
    // Get the track object and the track id
    const auto& track = inputTracks.getTrack(tid);

    // KF predicted momentum at the IP
    double me = 0.511 * Acts::UnitConstants::MeV;
    Acts::Vector3 pVec = track.momentum();
    double pMag = pVec.norm();
    m_ipMomentum.SetPxPyPzE(pVec.x(), pVec.y(), pVec.z(), std::hypot(pMag, me));

    // KF predicted IP momentum error
    m_ipMomentumError =
        TVector3(std::sqrt(track.covariance().diagonal().head<4>()[2]),
                 std::sqrt(track.covariance().diagonal().head<4>()[3]), 0);

    // KF predicted vertex position
    Acts::Vector3 vertex = {track.loc0(), 0, -track.loc1()};
    m_vertex = TVector3(vertex.x(), vertex.y(), vertex.z());

    // KF predicted vertex error
    Acts::Vector3 vertexError = {
        std::sqrt(track.covariance().diagonal().head<2>()[0]), 0,
        std::sqrt(track.covariance().diagonal().head<2>()[1])};
    m_vertexError = TVector3(vertexError.x(), vertexError.y(), vertexError.z());

    // Track hits from the measurements
    std::vector<TVector3> trackHits;

    // KF predicted track hits
    std::vector<TVector3> predictedTrackHits;
    std::vector<TVector3> filteredTrackHits;
    std::vector<TVector3> smoothedTrackHits;

    // KF residuals with respect to the measurements
    std::vector<TVector3> predictedResiduals;
    std::vector<TVector3> filteredResiduals;
    std::vector<TVector3> smoothedResiduals;

    // KF pulls with respect to the measurements
    std::vector<TVector3> predictedPulls;
    std::vector<TVector3> filteredPulls;
    std::vector<TVector3> smoothedPulls;

    double chi2Predicted = 0;
    double chi2Filtered = 0;
    double chi2Smoothed = 0;

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
      auto smoothedHit = state.effectiveProjector() * state.smoothed();

      auto hitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, hit, Acts::Vector3(1, 0, 0));
      auto predictedHitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, predictedHit, Acts::Vector3(1, 0, 0));
      auto filteredHitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, filteredHit, Acts::Vector3(1, 0, 0));
      auto smoothedHitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, smoothedHit, Acts::Vector3(1, 0, 0));

      // Get the residuals between the measurements and the predicted hits
      auto predictedResidual = hitGlobal - predictedHitGlobal;
      auto filteredResidual = hitGlobal - filteredHitGlobal;
      auto smoothedResidual = hitGlobal - smoothedHitGlobal;

      // KF predicted covariances
      auto measurementCov = state.effectiveCalibratedCovariance();

      // With respect to measurement
      auto predictedCov = state.effectiveProjector() *
                          state.predictedCovariance() *
                          state.effectiveProjector().transpose();

      auto filteredCov = state.effectiveProjector() *
                             state.filteredCovariance() *
                             state.effectiveProjector().transpose() -
                         measurementCov;

      auto smoothedCov = state.effectiveProjector() *
                             state.smoothedCovariance() *
                             state.effectiveProjector().transpose() -
                         measurementCov;

      auto predictedDiag =
          predictedCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
      auto filteredDiag =
          filteredCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
      auto smoothedDiag =
          smoothedCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

      // KF pulls with respect to the measurements
      auto predictedPull = predictedDiag.cwiseProduct(hit - predictedHit);
      auto filteredPull = filteredDiag.cwiseProduct(hit - filteredHit);
      auto smoothedPull = smoothedDiag.cwiseProduct(hit - smoothedHit);

      // Store the measurements hits
      trackHits.push_back(
          TVector3(hitGlobal.x(), hitGlobal.y(), hitGlobal.z()));

      // Store the KF predicted hits
      predictedTrackHits.push_back(TVector3(predictedHitGlobal.x(),
                                            predictedHitGlobal.y(),
                                            predictedHitGlobal.z()));
      filteredTrackHits.push_back(TVector3(
          filteredHitGlobal.x(), filteredHitGlobal.y(), filteredHitGlobal.z()));
      smoothedTrackHits.push_back(TVector3(
          smoothedHitGlobal.x(), smoothedHitGlobal.y(), smoothedHitGlobal.z()));

      // Store the residuals with respect to the measurements
      predictedResiduals.push_back(TVector3(
          predictedResidual.x(), predictedResidual.y(), predictedResidual.z()));
      filteredResiduals.push_back(TVector3(
          filteredResidual.x(), filteredResidual.y(), filteredResidual.z()));
      smoothedResiduals.push_back(TVector3(
          smoothedResidual.x(), smoothedResidual.y(), smoothedResidual.z()));

      // Store the pulls with respect to the measurements
      predictedPulls.push_back(
          TVector3(predictedPull.x(), 0, -predictedPull.y()));
      filteredPulls.push_back(TVector3(filteredPull.x(), 0, -filteredPull.y()));
      smoothedPulls.push_back(TVector3(smoothedPull.x(), 0, -smoothedPull.y()));

      // Add to the track chi2
      chi2Predicted += predictedPull.dot(predictedPull);
      chi2Filtered += filteredPull.dot(filteredPull);
      chi2Smoothed += smoothedPull.dot(smoothedPull);
    }

    // Measurement hits
    m_trackHits = trackHits;

    // KF predicted track hits
    m_predictedTrackHits = predictedTrackHits;
    m_filteredTrackHits = filteredTrackHits;
    m_smoothedTrackHits = smoothedTrackHits;

    // KF residuals with respect to the measurements
    m_predictedResiduals = predictedResiduals;
    m_filteredResiduals = filteredResiduals;
    m_smoothedResiduals = smoothedResiduals;

    // KF pulls with respect to the measurements
    m_predictedPulls = predictedPulls;
    m_filteredPulls = filteredPulls;
    m_smoothedPulls = smoothedPulls;

    // Chi2 of the tracks from different 
    // KF states
    m_chi2Predicted = chi2Predicted;
    m_chi2Filtered = chi2Filtered;
    m_chi2Smoothed = chi2Smoothed;

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
