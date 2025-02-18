#include "TrackingPipeline/Io/RootSimTrackCandidateWriter.hpp"

#include <algorithm>
#include <ranges>

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
  int buf_size = 32000;
  int split_lvl = 0;

  // True hits
  m_tree->Branch("trueTrackHits", &m_trueTrackHits, buf_size, split_lvl);

  // Measurement hits
  m_tree->Branch("trackHits", &m_trackHits, buf_size, split_lvl);

  // CKF predicted track hits
  m_tree->Branch("predictedTrackHits", &m_predictedTrackHits, buf_size,
                 split_lvl);
  m_tree->Branch("filteredTrackHits", &m_filteredTrackHits, buf_size,
                 split_lvl);

  // CKF residuals with respect to the true hits
  m_tree->Branch("truePredictedResiduals", &m_truePredictedResiduals, buf_size,
                 split_lvl);
  m_tree->Branch("trueFilteredResiduals", &m_trueFilteredResiduals, buf_size,
                 split_lvl);

  // CKF residuals with respect to the measurements
  m_tree->Branch("predictedResiduals", &m_predictedResiduals, buf_size,
                 split_lvl);
  m_tree->Branch("filteredResiduals", &m_filteredResiduals, buf_size,
                 split_lvl);

  // CKF pulls with respect to the true hits
  m_tree->Branch("truePredictedPulls", &m_truePredictedPulls, buf_size,
                 split_lvl);
  m_tree->Branch("trueFilteredPulls", &m_trueFilteredPulls, buf_size,
                 split_lvl);

  // CKF pulls with respect to the measurements
  m_tree->Branch("predictedPulls", &m_predictedPulls, buf_size, split_lvl);
  m_tree->Branch("filteredPulls", &m_filteredPulls, buf_size, split_lvl);

  // True momentum at the IP
  m_tree->Branch("ipMomentumTruth", &m_ipMomentumTruth);
  m_tree->Branch("vertexTruth", &m_vertexTruth);

  // Chi2 and ndf of the fitted track
  m_tree->Branch("chi2Predicted", &m_chi2Predicted, buf_size, split_lvl);
  m_tree->Branch("chi2Filtered", &m_chi2Filtered, buf_size, split_lvl);
  m_tree->Branch("ndf", &m_ndf, "ndf/I");

  // Matching degree between the true and the fitted track
  m_tree->Branch("matchingDegree", &m_matchingDegree, "matchingDegree/D");

  // Track ID
  m_tree->Branch("trackId", &m_trackId, "trackId/I");

  // Event ID
  m_tree->Branch("eventId", &m_eventId, "eventId/I");

  //------------------------------------------------------------------
  // Initialize the data handles
  m_inputTrackCandidates.initialize(m_cfg.inputTrackCandidates);
  m_inputTruthClusters.initialize(m_cfg.inputTruthClusters);
}

ProcessCode RootSimTrackCandidateWriter::finalize() {
  if (m_file) {
    m_file->Write();
    m_file->Close();
    delete m_file;
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootSimTrackCandidateWriter::write(const AlgorithmContext& ctx) {
  auto inputCandidates = m_inputTrackCandidates(ctx);

  auto inputTruthClusters = m_inputTruthClusters(ctx);

  std::lock_guard<std::mutex> lock(m_mutex);

  // Collect true track statistics
  std::map<TrackID, std::int32_t> trueTracksSig;

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

  m_truthSig = trueTracksSig.size();
  m_eventId = ctx.eventNumber;

  // Iterate over the fitted tracks
  for (std::size_t tid = 0; tid < inputCandidates.size(); tid++) {
    // Get the track object and the track id
    const auto& track = inputCandidates.getTrack(tid);

    // Track hits from the true information
    std::vector<TVector3> trueTrackHits;

    // Track hits from the measurements
    std::vector<TVector3> trackHits;

    // CKF predicted track hits
    std::vector<TVector3> predictedTrackHits;
    std::vector<TVector3> filteredTrackHits;

    // CKF residuals with respect to the true hits
    std::vector<TVector3> truePredictedResiduals;
    std::vector<TVector3> trueFilteredResiduals;

    // CKF residuals with respect to the measurements
    std::vector<TVector3> predictedResiduals;
    std::vector<TVector3> filteredResiduals;

    // CKF pulls with respect to the true hits
    std::vector<TVector3> truePredictedPulls;
    std::vector<TVector3> trueFilteredPulls;

    // CKF pulls with respect to the measurements
    std::vector<TVector3> predictedPulls;
    std::vector<TVector3> filteredPulls;

    // Flag indicating how many hits are matched
    // between the true and the fitted track
    double matchingDegree = 0;

    // Chi2 of different track states
    std::vector<double> chi2Predicted;
    std::vector<double> chi2Filtered;

    // Iterate over the track states
    std::map<TrackID, std::vector<std::int32_t>> trackStateIds;
    for (const auto& state : track.trackStatesReversed()) {
      // Skip the states without meaningful information
      if (!state.hasProjector()) {
        continue;
      }

      // Get the measurements source link
      auto sl = state.getUncalibratedSourceLink();
      auto ssl = sl.get<SimpleSourceLink>();

      auto cluster = inputTruthClusters.at(ssl.index());

      // Get the true hit
      Acts::Vector2 trueHit;
      TrackID currentTrackId;
      if (cluster.truthHits.size() == 0 || !cluster.isSignal) {
        trueHit = ssl.parameters();

        currentTrackId = std::make_tuple(-1, -1, -1);
      } else {
        auto sig = std::ranges::find_if(cluster.truthHits, [](const auto& hit) {
          return (hit.trackId == 1);
        });
        trueHit = sig->truthParameters.head<2>();

        currentTrackId =
            std::make_tuple(sig->trackId, sig->parentTrackId, sig->runId);
      }
      if (!trackStateIds.contains(currentTrackId)) {
        trackStateIds[currentTrackId] = {ssl.index()};
      } else {
        trackStateIds.at(currentTrackId).push_back(ssl.index());
      }

      // Get the true source link
      auto trueSl = Acts::SourceLink(cluster.sourceLink);

      // Get the measurements hit
      auto hit = state.effectiveCalibrated();

      // Project onto the prediction space
      auto predictedHit = state.effectiveProjector() * state.predicted();
      auto filteredHit = state.effectiveProjector() * state.filtered();

      // Transform the hits to the global coordinates
      auto trueHitGlobal = m_cfg.surfaceAccessor(trueSl)->localToGlobal(
          ctx.geoContext, trueHit, Acts::Vector3(1, 0, 0));
      auto hitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, hit, Acts::Vector3(1, 0, 0));
      auto predictedHitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, predictedHit, Acts::Vector3(1, 0, 0));
      auto filteredHitGlobal = state.referenceSurface().localToGlobal(
          ctx.geoContext, filteredHit, Acts::Vector3(1, 0, 0));

      // Get the residuals between the true and the predicted hits
      auto truePredictedResidual = trueHitGlobal - predictedHitGlobal;
      auto trueFilteredResidual = trueHitGlobal - filteredHitGlobal;

      // Get the residuals between the measurements and the predicted hits
      auto predictedResidual = hitGlobal - predictedHitGlobal;
      auto filteredResidual = hitGlobal - filteredHitGlobal;

      // CKF predicted covariances
      auto measurementCov = state.effectiveCalibratedCovariance();

      // With respect to truth
      auto predictedCovTruth =
          measurementCov + state.effectiveProjector() *
                               state.predictedCovariance() *
                               state.effectiveProjector().transpose();

      auto filteredCovTruth = state.effectiveProjector() *
                              state.filteredCovariance() *
                              state.effectiveProjector().transpose();

      // With respect to measurement
      auto predictedCov = state.effectiveProjector() *
                          state.predictedCovariance() *
                          state.effectiveProjector().transpose();

      auto filteredCov = state.effectiveProjector() *
                             state.filteredCovariance() *
                             state.effectiveProjector().transpose() -
                         measurementCov;

      // Extract diagonals
      auto predictedDiagTruth =
          predictedCovTruth.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
      auto filteredDiagTruth =
          filteredCovTruth.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

      auto predictedDiag =
          predictedCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
      auto filteredDiag =
          filteredCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

      // CKF pulls with respect to the true hits
      auto truePredictedPull =
          predictedDiagTruth.cwiseProduct(trueHit - predictedHit);
      auto trueFilteredPull =
          filteredDiagTruth.cwiseProduct(trueHit - filteredHit);

      // CKF pulls with respect to the measurements
      auto predictedPull = predictedDiag.cwiseProduct(hit - predictedHit);
      auto filteredPull = filteredDiag.cwiseProduct(hit - filteredHit);

      // Store the true hits
      trueTrackHits.push_back(
          TVector3(trueHitGlobal.x(), trueHitGlobal.y(), trueHitGlobal.z()));

      // Store the measurements hits
      trackHits.push_back(
          TVector3(hitGlobal.x(), hitGlobal.y(), hitGlobal.z()));

      // Store the KF predicted hits
      predictedTrackHits.push_back(TVector3(predictedHitGlobal.x(),
                                            predictedHitGlobal.y(),
                                            predictedHitGlobal.z()));
      filteredTrackHits.push_back(TVector3(
          filteredHitGlobal.x(), filteredHitGlobal.y(), filteredHitGlobal.z()));

      // Store the residuals with respect to the true hits
      truePredictedResiduals.push_back(TVector3(truePredictedResidual.x(),
                                                truePredictedResidual.y(),
                                                truePredictedResidual.z()));
      trueFilteredResiduals.push_back(TVector3(trueFilteredResidual.x(),
                                               trueFilteredResidual.y(),
                                               trueFilteredResidual.z()));

      // Store the residuals with respect to the measurements
      predictedResiduals.push_back(TVector3(
          predictedResidual.x(), predictedResidual.y(), predictedResidual.z()));
      filteredResiduals.push_back(TVector3(
          filteredResidual.x(), filteredResidual.y(), filteredResidual.z()));

      // Store the pulls with respect to the true hits
      truePredictedPulls.push_back(
          TVector3(truePredictedPull.x(), 0, -truePredictedPull.y()));
      trueFilteredPulls.push_back(
          TVector3(trueFilteredPull.x(), 0, -trueFilteredPull.y()));

      // Store the pulls with respect to the measurements
      predictedPulls.push_back(
          TVector3(predictedPull.x(), 0, -predictedPull.y()));
      filteredPulls.push_back(TVector3(filteredPull.x(), 0, -filteredPull.y()));

      // Chi2 of the track state
      chi2Predicted.push_back(predictedPull.dot(predictedPull));
      chi2Filtered.push_back(filteredPull.dot(filteredPull));
    }
    double me = 0.511 * Acts::UnitConstants::MeV;

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
    } else {
      // Get the true IP parameters
      std::int32_t refIndex = refTrackId->second.at(0);
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
          std::hypot(pivotHit->ipParameters.absoluteMomentum(), me));

      // Compute matching degree
      double trueTrackSize =
          std::ranges::count(trueTrackIds, refTrackId->first);

      if (trueTrackSize != m_cfg.targetTrueTrackSize) {
        continue;
      }

      matchingDegree = refTrackId->second.size() / trueTrackSize;
    }

    // True hits
    m_trueTrackHits = trueTrackHits;

    // Measurement hits
    m_trackHits = trackHits;

    // CKF predicted track hits
    m_predictedTrackHits = predictedTrackHits;
    m_filteredTrackHits = filteredTrackHits;

    // CKF residuals with respect to the true hits
    m_truePredictedResiduals = truePredictedResiduals;
    m_trueFilteredResiduals = trueFilteredResiduals;

    // CKF residuals with respect to the measurements
    m_predictedResiduals = predictedResiduals;
    m_filteredResiduals = filteredResiduals;

    // CKF pulls with respect to the true hits
    m_truePredictedPulls = truePredictedPulls;
    m_trueFilteredPulls = trueFilteredPulls;

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

    // Matching degree
    m_matchingDegree = matchingDegree;

    // Fill the tree
    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
