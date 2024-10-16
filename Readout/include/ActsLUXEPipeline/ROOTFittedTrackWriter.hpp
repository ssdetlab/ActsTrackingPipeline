#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"

#include "ActsLUXEPipeline/IWriter.hpp"
#include "ActsLUXEPipeline/ProcessCode.hpp"
#include "ActsLUXEPipeline/AlgorithmContext.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMatrix.h" 
#include "TLorentzVector.h" 

using namespace Acts::UnitLiterals;

using TrackID = std::tuple<std::int32_t, std::int32_t, std::int32_t>;

/// @brief Intermediate generalization of the 
/// ROOT file reader to be inhereted from by the
/// readers for the specific tree structures,
/// data types and geometries
///
/// @tparam measurementContainer_t container type 
/// for the measurements to be implemented 
///
/// @note The events are assumed to be ordered
// template <
// typename trajectory_t = Acts::VectorMultiTrajectory,
// typename container_t = Acts::VectorTrackContainer>
class ROOTFittedTrackWriter : public IWriter {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// Surface accessor
            Acts::SourceLinkSurfaceAccessor surfaceAccessor;
            /// Fitted track collection
            std::string inputKFTracks;
            /// Truth cluster data
            std::string inputTruthClusters;
            /// Name of the input tree
            std::string treeName;
            /// The names of the input files
            std::string filePath;
        };

        ROOTFittedTrackWriter(const ROOTFittedTrackWriter &) = delete;
        ROOTFittedTrackWriter(const ROOTFittedTrackWriter &&) = delete;
    
        /// Constructor
        /// @param config The Configuration struct
        ROOTFittedTrackWriter(const Config &config, Acts::Logging::Level level)
            : m_cfg(config),
            m_logger(Acts::getDefaultLogger(name(), level)) {
                if (m_cfg.filePath.empty()) {
                    throw std::invalid_argument("Missing filename");
                }
                if (m_cfg.treeName.empty()) {
                    throw std::invalid_argument("Missing tree name");
                }

                m_file = new TFile(m_cfg.filePath.c_str(), "RECREATE");
                m_trackTree = new TTree(m_cfg.treeName.c_str(), 
                    m_cfg.treeName.c_str());

                std::string flowTreeName = m_cfg.treeName + "_flow";
                m_flowTree = new TTree(flowTreeName.c_str(), 
                    flowTreeName.c_str());

                //------------------------------------------------------------------
                // Track tree branches
                int buf_size  = 32000;
                int split_lvl = 0;

                // True hits
                m_trackTree->Branch("trueTrackHits", &m_trueTrackHits, buf_size, split_lvl);
            
                // Measurement hits 
                m_trackTree->Branch("trackHits", &m_trackHits, buf_size, split_lvl);

                // KF predicted track hits
                m_trackTree->Branch("predictedTrackHits", &m_predictedTrackHits, buf_size, split_lvl);
                m_trackTree->Branch("filteredTrackHits", &m_filteredTrackHits, buf_size, split_lvl);
                m_trackTree->Branch("smoothedTrackHits", &m_smoothedTrackHits, buf_size, split_lvl);

                // KF residuals with respect to the true hits
                m_trackTree->Branch("truePredictedResiduals", &m_truePredictedResiduals, buf_size, split_lvl);
                m_trackTree->Branch("trueFilteredResiduals", &m_trueFilteredResiduals, buf_size, split_lvl);
                m_trackTree->Branch("trueSmoothedResiduals", &m_trueSmoothedResiduals, buf_size, split_lvl);

                // KF residuals with respect to the measurements
                m_trackTree->Branch("predictedResiduals", &m_predictedResiduals, buf_size, split_lvl);
                m_trackTree->Branch("filteredResiduals", &m_filteredResiduals, buf_size, split_lvl);
                m_trackTree->Branch("smoothedResiduals", &m_smoothedResiduals, buf_size, split_lvl);

                // KF pulls with respect to the true hits
                m_trackTree->Branch("truePredictedPulls", &m_truePredictedPulls, buf_size, split_lvl);
                m_trackTree->Branch("trueFilteredPulls", &m_trueFilteredPulls, buf_size, split_lvl);
                m_trackTree->Branch("trueSmoothedPulls", &m_trueSmoothedPulls, buf_size, split_lvl);

                // KF pulls with respect to the measurements
                m_trackTree->Branch("predictedPulls", &m_predictedPulls, buf_size, split_lvl);
                m_trackTree->Branch("filteredPulls", &m_filteredPulls, buf_size, split_lvl);
                m_trackTree->Branch("smoothedPulls", &m_smoothedPulls, buf_size, split_lvl);

                // KF predicted momentum at the IP
                m_trackTree->Branch("ipMomentum", &m_ipMomentum);
                m_trackTree->Branch("ipMomentumError", &m_ipMomentumError);
                m_trackTree->Branch("vertex", &m_vertex);
                m_trackTree->Branch("vertexError", &m_vertexError);

                // True momentum at the IP
                m_trackTree->Branch("ipMomentumTruth", &m_ipMomentumTruth);
                m_trackTree->Branch("vertexTruth", &m_vertexTruth);

                // Chi2 and ndf of the fitted track
                m_trackTree->Branch("chi2", &m_chi2, "chi2/D");
                m_trackTree->Branch("ndf", &m_ndf, "ndf/I");

                // Matching degree between the true and the fitted track
                m_trackTree->Branch("matchingDegree", &m_matchingDegree, "matchingDegree/D");

                // Track ID
                m_trackTree->Branch("trackId", &m_trackId, "trackId/I");

                // Event ID
                m_trackTree->Branch("eventId", &m_eventId, "eventId/I");

                //------------------------------------------------------------------
                // Flow tree branches
                m_flowTree->Branch("eventId", &m_eventId, "eventId/I");
                m_flowTree->Branch("TruthSig", &m_truthSig, "TruthSig/I");
                m_flowTree->Branch("KFReco", &m_KFReco, buf_size, split_lvl);

                //------------------------------------------------------------------
                // Initialize the data handles
                m_KFTracks.initialize(m_cfg.inputKFTracks);
                m_truthClusters.initialize(m_cfg.inputTruthClusters);
        }

        /// Destructor
        ~ROOTFittedTrackWriter() override {
            if (m_file) {
                m_file->Write();
                m_file->Close();
                delete m_file;
            }
        }
    
        /// Writer name() method
        std::string name() const { return "ROOTFittedTrackWriter"; }
    
        /// Write out data to the input stream
        ProcessCode write(const AlgorithmContext &ctx) override {
            auto inputKFTracks = m_KFTracks(ctx);

            auto inputTruthClusters = m_truthClusters(ctx);

            std::lock_guard<std::mutex> lock(m_mutex);

            // Collect true track statistics
            std::map<TrackID, std::int32_t> trueTracksSig;
            std::map<TrackID, std::int32_t> trueTracksBkg;
            for (auto& cluster : inputTruthClusters) {
                for (int i = 0; i < cluster.truthHits.size(); i++) {
                    SimHit truthHit = cluster.truthHits.at(i);

                    TrackID trueTrackId = 
                        std::make_tuple(
                            truthHit.trackId, 
                            truthHit.parentTrackId,
                            truthHit.runId);
                    
                    if (truthHit.trackId == 1) {
                        if (trueTracksSig.find(trueTrackId) == trueTracksSig.end()) {
                            trueTracksSig.insert({trueTrackId, 1});
                        } else {
                            trueTracksSig.at(trueTrackId)++;
                        }
                    }
                    else {
                        if (trueTracksBkg.find(trueTrackId) == trueTracksBkg.end()) {
                            trueTracksBkg.insert({trueTrackId, 1});
                        } else {
                            trueTracksBkg.at(trueTrackId)++;
                        }
                    }
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

                // Iterate over the track states
                std::vector<std::pair<TrackID, std::int32_t>> trackStateIds;
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
                    for (int i = 0; i < cluster.truthHits.size(); i++) {
                        SimHit truthHit = cluster.truthHits.at(i);
                        if (!cluster.isSignal || truthHit.trackId == 1) {
                            trueHit = truthHit.truthParameters.head<2>();
                            
                            TrackID currentTrackId = 
                                std::make_tuple(
                                    truthHit.trackId, 
                                    truthHit.parentTrackId,
                                    truthHit.runId);

                            trackStateIds.push_back(
                                {currentTrackId, ssl.index()});
                            break;
                        }
                    }
                    if (cluster.truthHits.size() == 0) {
                        trueHit = ssl.parameters();

                        TrackID currentTrackId = 
                            std::make_tuple(
                                -1, 
                                -1,
                                -1);

                        trackStateIds.push_back(
                            {currentTrackId, ssl.index()});
                    }

                    // Get the true source link
                    auto trueSl = Acts::SourceLink(cluster.sourceLink);

                    // Get the measurements hit
                    auto hit = state.effectiveCalibrated();

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
                    auto predictedCovariance =
                        state.effectiveProjector() * 
                        state.predictedCovariance() * 
                        state.effectiveProjector().transpose();

                    auto filteredCovariance =
                        state.effectiveProjector() * 
                        state.filteredCovariance() * 
                        state.effectiveProjector().transpose();

                    auto smoothedCovariance =
                        state.effectiveProjector() * 
                        state.smoothedCovariance() * 
                        state.effectiveProjector().transpose();

                    // KF pulls with respect to the true hits
                    auto truePredictedPull = predictedCovariance.inverse().cwiseSqrt() * (trueHit - predictedHit);
                    auto trueFilteredPull = filteredCovariance.inverse().cwiseSqrt() * (trueHit - filteredHit);
                    auto trueSmoothedPull = smoothedCovariance.inverse().cwiseSqrt() * (trueHit - smoothedHit);

                    // KF pulls with respect to the measurements
                    auto predictedPull = predictedCovariance.inverse().cwiseSqrt() * (hit - predictedHit);
                    auto filteredPull = filteredCovariance.inverse().cwiseSqrt() * (hit - filteredHit);
                    auto smoothedPull = smoothedCovariance.inverse().cwiseSqrt() * (hit - smoothedHit);

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
                        TVector3(predictedPull.x(), 0, -predictedPull.y()));
                    trueFilteredPulls.push_back(
                        TVector3(filteredPull.x(), 0, -filteredPull.y()));
                    trueSmoothedPulls.push_back(
                        TVector3(smoothedPull.x(), 0, -smoothedPull.y()));

                    // Store the pulls with respect to the measurements
                    predictedPulls.push_back(
                        TVector3(predictedPull.x(), 0, -predictedPull.y()));
                    filteredPulls.push_back(
                        TVector3(filteredPull.x(), 0, -filteredPull.y()));
                    smoothedPulls.push_back(
                        TVector3(smoothedPull.x(), 0, -smoothedPull.y()));
                }
                if (trackStateIds.size() != 4) {
                    throw std::runtime_error("Track does not have 4 hits");
                }

                // Count the number of track ID occurrences
                std::map<TrackID, std::int32_t> trackIDcount;
                for (auto [currentTrackId, index] : trackStateIds) {
                    if (trackIDcount.find(currentTrackId) == trackIDcount.end()) {
                        trackIDcount.insert({currentTrackId, 1});
                    } else {
                        trackIDcount.at(currentTrackId)++;
                    }
                }
                // Find the track ID with the highest occurrence rate
                TrackID refTrackId;
                std::int32_t refCount = 0;
                for (auto [currentTrackId, count] : trackIDcount) {
                    if (count > refCount && std::get<0>(currentTrackId) == 1) {
                        refTrackId = currentTrackId;
                        refCount = count;
                    }
                }
                if (refCount == 0) {
                    matchingDegree = 0;
                }
                else {
                    // Extract the index of the reference track ID
                    std::int32_t refIndex = 0;
                    for (auto [currentTrackId, index] : trackStateIds) {
                        if (currentTrackId == refTrackId) {
                            refIndex = index;
                            break;
                        }
                    }
    
                    // Get the true IP parameters
                    auto cluster = inputTruthClusters.at(refIndex);
                    for (int j = 0; j < cluster.truthHits.size(); j++) {
                        SimHit truthHit = cluster.truthHits.at(j);
    
                        TrackID pivotTrackId = std::make_tuple(
                            truthHit.trackId, 
                            truthHit.parentTrackId,
                            truthHit.runId);
    
                        if (pivotTrackId == refTrackId) {
                            m_ipMomentumTruth.SetPxPyPzE(
                                truthHit.ipParameters.momentum().x(),
                                truthHit.ipParameters.momentum().y(),
                                truthHit.ipParameters.momentum().z(),
                                std::hypot(truthHit.ipParameters.absoluteMomentum(), me));
                        }
                    }
                    for (auto [currentTrackId, index] : trackStateIds) {
                        if (currentTrackId == refTrackId) {
                            matchingDegree += 1;
                        }
                    }
                    matchingDegree /= trackStateIds.size();
                }

                if (m_KFReco.find(matchingDegree) == m_KFReco.end()) {
                    m_KFReco.insert({matchingDegree, 1});
                } else {
                    m_KFReco.at(matchingDegree)++;
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
                m_chi2 = track.chi2();

                // Number of degrees of freedom
                m_ndf = track.nDoF();

                // Track Id
                m_trackId = id;

                // Matching degree
                m_matchingDegree = matchingDegree;

                // Fill the tree
                m_trackTree->Fill();
            }

            m_flowTree->Fill();

            m_KFReco.clear();

            // Return success flag
            return ProcessCode::SUCCESS;
        }
    
        /// Readonly access to the config
        const Config &config() const { return m_cfg; }

    private:
        /// Private access to the logging instance
        const Acts::Logger &logger() const { return *m_logger; }

        /// The config class
        Config m_cfg;

        ReadDataHandle<Tracks<
            Acts::VectorTrackContainer,
            Acts::VectorMultiTrajectory>>
                m_KFTracks{this, "KFTracks"};  

        ReadDataHandle<SimClusters> m_truthClusters{this, "TruthClusters"};


        std::unique_ptr<const Acts::Logger> m_logger;

        /// The output file
        TFile *m_file = nullptr;

        /// The output tree
        TTree *m_trackTree = nullptr;

        /// Cut flow tree
        TTree *m_flowTree = nullptr;

    protected:
        /// True hits
        std::vector<TVector3> m_trueTrackHits;

        /// Measurement hits
        std::vector<TVector3> m_trackHits;
        
        /// KF predicted track hits
        std::vector<TVector3> m_predictedTrackHits;
        std::vector<TVector3> m_filteredTrackHits;
        std::vector<TVector3> m_smoothedTrackHits;

        /// KF residuals with respect to the true hits
        std::vector<TVector3> m_truePredictedResiduals;
        std::vector<TVector3> m_trueFilteredResiduals;
        std::vector<TVector3> m_trueSmoothedResiduals;

        /// KF residuals with respect to the measurements
        std::vector<TVector3> m_predictedResiduals;
        std::vector<TVector3> m_filteredResiduals;
        std::vector<TVector3> m_smoothedResiduals;

        /// KF pulls with respect to the true hits
        std::vector<TVector3> m_truePredictedPulls;
        std::vector<TVector3> m_trueFilteredPulls;
        std::vector<TVector3> m_trueSmoothedPulls;

        /// KF pulls with respect to the measurements
        std::vector<TVector3> m_predictedPulls;
        std::vector<TVector3> m_filteredPulls;
        std::vector<TVector3> m_smoothedPulls;

        /// Chi2 of the track
        double m_chi2;

        /// Number of degrees of freedom
        /// of the track
        int m_ndf;

        /// Matching degree
        double m_matchingDegree;

        /// TrackId
        int m_trackId;

        /// EventId
        int m_eventId;

        /// KF predicted momentum at the IP
        TLorentzVector m_ipMomentum;
        TVector3 m_ipMomentumError;
        TVector3 m_vertex;
        TVector3 m_vertexError;

        /// True momentum at the IP
        TLorentzVector m_ipMomentumTruth;
        TVector3 m_vertexTruth;

        /// Number of true tracks prior to 
        /// applying the cuts
        std::int32_t m_truthSig;

        /// Number of KF reconstructed tracks
        /// sorted by the matching degree
        std::map<double, std::int32_t> m_KFReco;

        /// Mutex to protect the tree filling
        std::mutex m_mutex;
};
