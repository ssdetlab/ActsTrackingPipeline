#include "TrackingPipeline/Io/RootFittedSimTrackWriter.hpp"

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

#include <algorithm>
#include <ranges>

RootFittedSimTrackWriter::RootFittedSimTrackWriter(
    const Config& config, Acts::Logging::Level level)
        : m_cfg(config),
        m_logger(Acts::getDefaultLogger(name(), level)) {
            if (m_cfg.filePath.empty()) {
                throw std::invalid_argument("Missing filename");
            }
            if (m_cfg.treeName.empty()) {
                throw std::invalid_argument("Missing tree name");
            }

            m_file = new TFile(m_cfg.filePath.c_str(), "RECREATE");
            m_tree = new TTree(m_cfg.treeName.c_str(), 
                m_cfg.treeName.c_str());

            //------------------------------------------------------------------
            // Track tree branches
            int buf_size  = 32000;
            int split_lvl = 0;

            // True hits
            m_tree->Branch("trueTrackHits", &m_trueTrackHits, buf_size, split_lvl);
        
            // Measurement hits 
            m_tree->Branch("trackHits", &m_trackHits, buf_size, split_lvl);

            // KF predicted track hits
            m_tree->Branch("predictedTrackHits", &m_predictedTrackHits, buf_size, split_lvl);
            m_tree->Branch("filteredTrackHits", &m_filteredTrackHits, buf_size, split_lvl);
            m_tree->Branch("smoothedTrackHits", &m_smoothedTrackHits, buf_size, split_lvl);

            // KF residuals with respect to the true hits
            m_tree->Branch("truePredictedResiduals", &m_truePredictedResiduals, buf_size, split_lvl);
            m_tree->Branch("trueFilteredResiduals", &m_trueFilteredResiduals, buf_size, split_lvl);
            m_tree->Branch("trueSmoothedResiduals", &m_trueSmoothedResiduals, buf_size, split_lvl);

            // KF residuals with respect to the measurements
            m_tree->Branch("predictedResiduals", &m_predictedResiduals, buf_size, split_lvl);
            m_tree->Branch("filteredResiduals", &m_filteredResiduals, buf_size, split_lvl);
            m_tree->Branch("smoothedResiduals", &m_smoothedResiduals, buf_size, split_lvl);

            // KF pulls with respect to the true hits
            m_tree->Branch("truePredictedPulls", &m_truePredictedPulls, buf_size, split_lvl);
            m_tree->Branch("trueFilteredPulls", &m_trueFilteredPulls, buf_size, split_lvl);
            m_tree->Branch("trueSmoothedPulls", &m_trueSmoothedPulls, buf_size, split_lvl);

            // KF pulls with respect to the measurements
            m_tree->Branch("predictedPulls", &m_predictedPulls, buf_size, split_lvl);
            m_tree->Branch("filteredPulls", &m_filteredPulls, buf_size, split_lvl);
            m_tree->Branch("smoothedPulls", &m_smoothedPulls, buf_size, split_lvl);

            // KF predicted momentum at the IP
            m_tree->Branch("ipMomentum", &m_ipMomentum);
            m_tree->Branch("ipMomentumError", &m_ipMomentumError);
            m_tree->Branch("vertex", &m_vertex);
            m_tree->Branch("vertexError", &m_vertexError);

            // True momentum at the IP
            m_tree->Branch("ipMomentumTruth", &m_ipMomentumTruth);
            m_tree->Branch("vertexTruth", &m_vertexTruth);

            // Chi2 and ndf of the fitted track
            m_tree->Branch("chi2", &m_chi2, "chi2/D");
            m_tree->Branch("ndf", &m_ndf, "ndf/I");

            // Matching degree between the true and the fitted track
            m_tree->Branch("matchingDegree", &m_matchingDegree, "matchingDegree/D");

            // Track ID
            m_tree->Branch("trackId", &m_trackId, "trackId/I");

            // Event ID
            m_tree->Branch("eventId", &m_eventId, "eventId/I");

            //------------------------------------------------------------------
            // Initialize the data handles
            m_KFTracks.initialize(m_cfg.inputKFTracks);
            m_truthClusters.initialize(m_cfg.inputTruthClusters);
}

ProcessCode RootFittedSimTrackWriter::finalize() {
    if (m_file) {
        m_file->Write();
        m_file->Close();
        delete m_file;
    }
    return ProcessCode::SUCCESS;
}

ProcessCode RootFittedSimTrackWriter::write(const AlgorithmContext &ctx) {
    auto inputKFTracks = m_KFTracks(ctx);

    auto inputTruthClusters = m_truthClusters(ctx);

    std::lock_guard<std::mutex> lock(m_mutex);

    // Collect true track statistics
    std::map<TrackID, std::int32_t> trueTracksSig;

    // Collect true track statistics
    auto trueTrackIds = inputTruthClusters 
        | std::views::filter([](const auto& cl) {
                return cl.isSignal;
            })
        | std::views::transform([](const auto& cl) {
            return cl.truthHits;
        })
        | std::views::join
        | std::views::transform([](const auto& hit) 
            -> TrackID {
            return {hit.trackId, hit.parentTrackId, hit.runId};
        });
    
    for (const auto& id : trueTrackIds) {
        if (!trueTracksSig.contains(id)) {
            trueTracksSig[id] = 1;
        }
        else {
            trueTracksSig.at(id)++;
        }
    }

    m_truthSig = trueTracksSig.size();
    m_eventId = ctx.eventNumber;

    // Iterate over the fitted tracks
    for (int idx = 0; idx < inputKFTracks.size(); idx++) {
        // Get the track object and the track id
        auto [id,track] = inputKFTracks.getByIndex(idx);

        // KF predicted momentum at the IP
        double me = 0.511 * Acts::UnitConstants::MeV;
        Acts::Vector3 pVec = track.momentum();
        double pMag = pVec.norm();
        m_ipMomentum.SetPxPyPzE(
            pVec.x(), pVec.y(), pVec.z(), std::hypot(pMag, me));

        // KF predicted IP momentum error
        m_ipMomentumError = TVector3(
            std::sqrt(track.covariance().diagonal().head<4>()[2]),
            std::sqrt(track.covariance().diagonal().head<4>()[3]),
            0);

        // KF predicted vertex position
        Acts::Vector3 vertex = {track.loc0(), 0, -track.loc1()};
        m_vertex = TVector3(vertex.x(), vertex.y(), vertex.z());

        // KF predicted vertex error
        Acts::Vector3 vertexError = {
            std::sqrt(track.covariance().diagonal().head<2>()[0]),
            0,
            std::sqrt(track.covariance().diagonal().head<2>()[1])};
        m_vertexError = TVector3(vertexError.x(), vertexError.y(), vertexError.z());

        // Track hits from the true information
        std::vector<TVector3> trueTrackHits;

        // Track hits from the measurements
        std::vector<TVector3> trackHits;

        // KF predicted track hits
        std::vector<TVector3> predictedTrackHits;
        std::vector<TVector3> filteredTrackHits;
        std::vector<TVector3> smoothedTrackHits;

        // KF residuals with respect to the true hits
        std::vector<TVector3> truePredictedResiduals;
        std::vector<TVector3> trueFilteredResiduals;
        std::vector<TVector3> trueSmoothedResiduals;

        // KF residuals with respect to the measurements
        std::vector<TVector3> predictedResiduals;
        std::vector<TVector3> filteredResiduals;
        std::vector<TVector3> smoothedResiduals;

        // KF pulls with respect to the true hits
        std::vector<TVector3> truePredictedPulls;
        std::vector<TVector3> trueFilteredPulls;
        std::vector<TVector3> trueSmoothedPulls;

        // KF pulls with respect to the measurements
        std::vector<TVector3> predictedPulls;
        std::vector<TVector3> filteredPulls;
        std::vector<TVector3> smoothedPulls;
        
        // Flag indicating how many hits are matched
        // between the true and the fitted track
        double matchingDegree = 0;

        double trueChi2 = 0;

        // Iterate over the track states
        std::map<TrackID, std::vector<std::int32_t>> trackStateIds;
        for (auto state : track.trackStatesReversed()) {
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

                currentTrackId = 
                    std::make_tuple(
                        -1, 
                        -1,
                        -1);
            }
            else {
                auto sig = std::ranges::find_if(
                    cluster.truthHits, 
                    [](const auto& hit) {
                        return (hit.trackId == 1);
                });
                trueHit = sig->truthParameters.head<2>();

                currentTrackId = 
                    std::make_tuple(
                        sig->trackId, 
                        sig->parentTrackId,
                        sig->runId);
            }
            if (!trackStateIds.contains(currentTrackId)) {
                trackStateIds[currentTrackId] = {ssl.index()};
            }
            else {
                trackStateIds.at(currentTrackId).push_back(ssl.index());
            }

            // Get the true source link
            auto trueSl = Acts::SourceLink(cluster.sourceLink);

            // Get the measurements hit
            auto hit = state.effectiveCalibrated();

            // std::cout << "HIT: " << hit.transpose() << "\n";

            // Project onto the prediction space
            auto predictedHit = state.effectiveProjector() * state.predicted();
            auto filteredHit = state.effectiveProjector() * state.filtered();
            auto smoothedHit = state.effectiveProjector() * state.smoothed();

            // Transform the hits to the global coordinates
            auto trueHitGlobal = m_cfg.surfaceAccessor(trueSl)->localToGlobal(
                ctx.geoContext, trueHit, Acts::Vector3(1, 0, 0));
            auto hitGlobal = state.referenceSurface().localToGlobal(
                ctx.geoContext, hit, Acts::Vector3(1, 0, 0));
            auto predictedHitGlobal = state.referenceSurface().localToGlobal(
                ctx.geoContext, predictedHit, Acts::Vector3(1, 0, 0));
            auto filteredHitGlobal = state.referenceSurface().localToGlobal(
                ctx.geoContext, filteredHit, Acts::Vector3(1, 0, 0));
            auto smoothedHitGlobal = state.referenceSurface().localToGlobal(
                ctx.geoContext, smoothedHit, Acts::Vector3(1, 0, 0));

            // Get the residuals between the true and the predicted hits
            auto truePredictedResidual = trueHitGlobal - predictedHitGlobal;
            auto trueFilteredResidual = trueHitGlobal - filteredHitGlobal;
            auto trueSmoothedResidual = trueHitGlobal - smoothedHitGlobal;

            // Get the residuals between the measurements and the predicted hits
            auto predictedResidual = hitGlobal - predictedHitGlobal;
            auto filteredResidual = hitGlobal - filteredHitGlobal;
            auto smoothedResidual = hitGlobal - smoothedHitGlobal;

            // KF predicted covariances
            auto measurementCov = 
                state.effectiveCalibratedCovariance();

            // With respect to truth
            auto predictedCovTruth = measurementCov +  
                state.effectiveProjector() * 
                state.predictedCovariance() * 
                state.effectiveProjector().transpose();

            auto filteredCovTruth =
                state.effectiveProjector() * 
                state.filteredCovariance() * 
                state.effectiveProjector().transpose();

            auto smoothedCovTruth =
                state.effectiveProjector() * 
                state.smoothedCovariance() * 
                state.effectiveProjector().transpose();

            // With respect to measurement
            auto predictedCov =
                state.effectiveProjector() * 
                state.predictedCovariance() * 
                state.effectiveProjector().transpose();

            auto filteredCov =
                state.effectiveProjector() * 
                state.filteredCovariance() * 
                state.effectiveProjector().transpose() - measurementCov;

            auto smoothedCov =
                state.effectiveProjector() * 
                state.smoothedCovariance() * 
                state.effectiveProjector().transpose() - measurementCov;

            // Extract diagonals
            auto predictedDiagTruth = 
                predictedCovTruth.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
            auto filteredDiagTruth = 
                filteredCovTruth.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
            auto smoothedDiagTruth = 
                smoothedCovTruth.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

            auto predictedDiag = 
                predictedCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
            auto filteredDiag = 
                filteredCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
            auto smoothedDiag = 
                smoothedCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

            // KF pulls with respect to the true hits
            auto truePredictedPull = predictedDiagTruth.cwiseProduct(trueHit - predictedHit);
            auto trueFilteredPull = filteredDiagTruth.cwiseProduct(trueHit - filteredHit);
            auto trueSmoothedPull = smoothedDiagTruth.cwiseProduct(trueHit - smoothedHit);

            // KF pulls with respect to the measurements
            auto predictedPull = predictedDiag.cwiseProduct(hit - predictedHit);
            auto filteredPull = filteredDiag.cwiseProduct(hit - filteredHit);
            auto smoothedPull = smoothedDiag.cwiseProduct(hit - smoothedHit);

            // Store the true hits
            trueTrackHits.push_back(
                TVector3(trueHitGlobal.x(), trueHitGlobal.y(), trueHitGlobal.z()));

            // Store the measurements hits
            trackHits.push_back(
                TVector3(hitGlobal.x(), hitGlobal.y(), hitGlobal.z()));

            // Store the KF predicted hits
            predictedTrackHits.push_back(
                TVector3(predictedHitGlobal.x(), predictedHitGlobal.y(), predictedHitGlobal.z()));
            filteredTrackHits.push_back(
                TVector3(filteredHitGlobal.x(), filteredHitGlobal.y(), filteredHitGlobal.z()));
            smoothedTrackHits.push_back(
                TVector3(smoothedHitGlobal.x(), smoothedHitGlobal.y(), smoothedHitGlobal.z()));

            // Store the residuals with respect to the true hits
            truePredictedResiduals.push_back(
                TVector3(truePredictedResidual.x(), truePredictedResidual.y(), truePredictedResidual.z()));
            trueFilteredResiduals.push_back(
                TVector3(trueFilteredResidual.x(), trueFilteredResidual.y(), trueFilteredResidual.z()));
            trueSmoothedResiduals.push_back(
                TVector3(trueSmoothedResidual.x(), trueSmoothedResidual.y(), trueSmoothedResidual.z()));

            // Store the residuals with respect to the measurements
            predictedResiduals.push_back(
                TVector3(predictedResidual.x(), predictedResidual.y(), predictedResidual.z()));
            filteredResiduals.push_back(
                TVector3(filteredResidual.x(), filteredResidual.y(), filteredResidual.z()));
            smoothedResiduals.push_back(
                TVector3(smoothedResidual.x(), smoothedResidual.y(), smoothedResidual.z()));

            // Store the pulls with respect to the true hits
            truePredictedPulls.push_back(
                TVector3(truePredictedPull.x(), 0, -truePredictedPull.y()));
            trueFilteredPulls.push_back(
                TVector3(trueFilteredPull.x(), 0, -trueFilteredPull.y()));
            trueSmoothedPulls.push_back(
                TVector3(trueSmoothedPull.x(), 0, -trueSmoothedPull.y()));

            // Store the pulls with respect to the measurements
            predictedPulls.push_back(
                TVector3(predictedPull.x(), 0, -predictedPull.y()));
            filteredPulls.push_back(
                TVector3(filteredPull.x(), 0, -filteredPull.y()));
            smoothedPulls.push_back(
                TVector3(smoothedPull.x(), 0, -smoothedPull.y()));
        }

        // Matching degree is computed with respect 
        // to the most often occuring signal track
        auto refTrackId = std::ranges::max_element(
            trackStateIds,
            [](const auto& pairA, const auto& pairB) {
                if (std::get<0>(pairA.first) == -1 
                    && std::get<0>(pairB.first) != -1) {
                        return true;
                }
                else if (std::get<0>(pairA.first) != -1 
                    && std::get<0>(pairB.first) == -1) {
                        return false;
                }
                else {
                    return pairA.second.size() < pairB.second.size();
                }

        });
        if (std::get<0>(refTrackId->first) == -1) {
            matchingDegree = 0;
        }
        else {
            // Get the true IP parameters
            std::int32_t refIndex = refTrackId->second.at(0);
            auto cluster = inputTruthClusters.at(refIndex);
            auto pivotHit = std::ranges::find_if(
                cluster.truthHits,
                [&](const auto& hit) {
                    TrackID id{hit.trackId, hit.parentTrackId, hit.runId};
                    return (id == refTrackId->first);
            });

            m_ipMomentumTruth.SetPxPyPzE(
                pivotHit->ipParameters.momentum().x(),
                pivotHit->ipParameters.momentum().y(),
                pivotHit->ipParameters.momentum().z(),
                std::hypot(pivotHit->ipParameters.absoluteMomentum(), me));

            // Compute matching degree
            double trueTrackSize = std::ranges::count(
                trueTrackIds, refTrackId->first);
            
            if (trueTrackSize != m_cfg.targetTrueTrackSize) {
                continue;
            }

            matchingDegree = refTrackId->second.size() / trueTrackSize;
        }

        // True hits
        m_trueTrackHits = trueTrackHits;

        // Measurement hits
        m_trackHits = trackHits;
        
        // KF predicted track hits
        m_predictedTrackHits = predictedTrackHits;
        m_filteredTrackHits = filteredTrackHits;
        m_smoothedTrackHits = smoothedTrackHits;
        
        // KF residuals with respect to the true hits
        m_truePredictedResiduals = truePredictedResiduals;
        m_trueFilteredResiduals = trueFilteredResiduals;
        m_trueSmoothedResiduals = trueSmoothedResiduals;

        // KF residuals with respect to the measurements
        m_predictedResiduals = predictedResiduals;
        m_filteredResiduals = filteredResiduals;
        m_smoothedResiduals = smoothedResiduals;

        // KF pulls with respect to the true hits
        m_truePredictedPulls = truePredictedPulls;
        m_trueFilteredPulls = trueFilteredPulls;
        m_trueSmoothedPulls = trueSmoothedPulls;

        // KF pulls with respect to the measurements
        m_predictedPulls = predictedPulls;
        m_filteredPulls = filteredPulls;
        m_smoothedPulls = smoothedPulls;

        // Chi2 of the track
        // with respect ot the 
        // measurement
        m_chi2 = track.chi2();

        // Number of degrees of freedom
        m_ndf = track.nDoF();

        // Track Id
        m_trackId = id;

        // Matching degree
        m_matchingDegree = matchingDegree;

        // Fill the tree
        m_tree->Fill();
    }

    // Return success flag
    return ProcessCode::SUCCESS;
}
