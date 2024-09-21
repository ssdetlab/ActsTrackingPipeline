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
            /// Name of the fitted track collection
            std::string inputTrackCollection;
            /// Name of the seed collection
            std::string inputSeedCollection;
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
                m_tree = new TTree(m_cfg.treeName.c_str(), 
                    m_cfg.treeName.c_str());

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

                // KF distances with respect to the true hits
                m_tree->Branch("truePredictedDistances", &m_truePredictedDistances, buf_size, split_lvl);
                m_tree->Branch("trueFilteredDistances", &m_trueFilteredDistances, buf_size, split_lvl);
                m_tree->Branch("trueSmoothedDistances", &m_trueSmoothedDistances, buf_size, split_lvl);

                // KF distances with respect to the measurements
                m_tree->Branch("predictedDistances", &m_predictedDistances, buf_size, split_lvl);
                m_tree->Branch("filteredDistances", &m_filteredDistances, buf_size, split_lvl);
                m_tree->Branch("smoothedDistances", &m_smoothedDistances, buf_size, split_lvl);

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

                // Flag indicating if the number of measurements
                // and the number of true source links are different
                m_tree->Branch("mismatch", &m_mismatch, "mismatch/I");

                // Track ID
                m_tree->Branch("trackId", &m_trackId, "trackId/I");

                m_inputTracks.initialize(m_cfg.inputTrackCollection);
                m_inputSeeds.initialize(m_cfg.inputSeedCollection);
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
            auto inputTracks = m_inputTracks(ctx);

            auto inputSeeds = m_inputSeeds(ctx);

            std::lock_guard<std::mutex> lock(m_mutex);

            int nSkippedStates = 0;

            // Iterate over the fitted tracks
            for (int idx = 0; idx < inputTracks.size(); idx++) {
                // Get the track object and the track id
                auto [id,track] = inputTracks.getByIndex(idx);

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

                // Iterate over the true seed information
                // and find the true source links
                std::vector<Acts::SourceLink> trueSourceLinks;
                for (auto& seed : inputSeeds) {
                    if (seed.trackId == id) {
                        for (auto& sl : seed.sourceLinks) {
                            trueSourceLinks.push_back(sl);
                        }

                        // True momentum at the IP
                        Acts::Vector3 pVecTruth = seed.ipParameters.momentum();
                        double pMagTruth = pVecTruth.norm();
                        m_ipMomentumTruth.SetPxPyPzE(
                            pVecTruth.x(), pVecTruth.y(), pVecTruth.z(), std::hypot(pMagTruth, me));

                        // True vertex position                        
                        Acts::Vector3 vertexTruth = {
                            seed.ipParameters.position().x(),
                            seed.ipParameters.position().y(),
                            seed.ipParameters.position().z()};
                        m_vertexTruth = TVector3(vertexTruth.x(), vertexTruth.y(), vertexTruth.z());

                        break;
                    }
                }

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

                // KF distances with respect to the true hits
                std::vector<double> truePredictedDistances;
                std::vector<double> trueFilteredDistances;
                std::vector<double> trueSmoothedDistances;

                // KF distances with respect to the measurements
                std::vector<double> predictedDistances;
                std::vector<double> filteredDistances;
                std::vector<double> smoothedDistances;

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

                // Flag indicating if the number of measurements
                // and the number of true source links are different
                m_mismatch = 0;

                // Order true source links by the order of the track states
                std::sort(trueSourceLinks.begin(), trueSourceLinks.end(),
                    [](const Acts::SourceLink& a, const Acts::SourceLink& b) {
                        std::uint64_t aLayerId = 
                            a.get<SimpleSourceLink>().geometryId().sensitive() / 10;
                        std::uint64_t bLayerId = 
                            b.get<SimpleSourceLink>().geometryId().sensitive() / 10;
                        return aLayerId < bLayerId;
                    });

                int norm = trueSourceLinks.size();

                // Iterate over the track states
                for (auto state : track.trackStatesReversed()) {
                    // Skip the states without meaningful information
                    if (!state.hasProjector()) {
                        continue;
                    }
                    if (trueSourceLinks.empty()) {
                        break;
                    }
                    // Get the true source link
                    auto trueSl = trueSourceLinks.back();
                    auto trueSsl = trueSl.get<SimpleSourceLink>();

                    // Get the measurements source link
                    auto sl = state.getUncalibratedSourceLink();
                    auto ssl = sl.get<SimpleSourceLink>();
                    
                    // Check if the source link is in the true source links
                    if (trueSsl == ssl) {
                        matchingDegree++;
                        trueSourceLinks.pop_back();
                    }
                    else if (track.nMeasurements() == norm) {
                        // Mismatch in the same layer
                        trueSourceLinks.pop_back();
                    }
                    else if (track.nMeasurements() > norm) {
                        // Virtually means that there are three true source links
                        // and four measurements passed to the KF
                        // True source link can either occupy the last three
                        // layers (mismatched the pivot) or the first three
                        // layers (fake hit in the last layer)
                        
                        std::uint64_t trueLayerId = 
                            trueSsl.geometryId().sensitive() / 10;
                        std::uint64_t layerId =
                            ssl.geometryId().sensitive() / 10;
                        
                        // If the layer id of the true source link is less than
                        // the layer id of the measurement source link
                        // it means that the last measurement is a fake hit
                        // Skip the last measurement
                        if (trueLayerId < layerId) {
                            m_mismatch = 1;
                            continue;
                        }
                        // If the layer id of the true source link is greater than
                        // the layer id of the measurement source link
                        // it means that the pivot is mismatched
                        else if (trueLayerId > layerId) {
                            trueSourceLinks.pop_back();
                            m_mismatch = 2;
                            continue;
                        }
                        // Do nothing if the mismatch happened in the same layer
                        // as it is a regular CKF mismatch
                        trueSourceLinks.pop_back();
                    }
                    else if (track.nMeasurements() < norm) {
                        // Virtually means that there are four true source links
                        // and three measurements passed to the KF
                        // Measurements can either occupy the first three
                        // layers (missed the last measurement) or the last three
                        // layers (mismatched the pivot)
                        
                        std::uint64_t trueLayerId = 
                            trueSsl.geometryId().sensitive() / 10;
                        std::uint64_t layerId =
                            ssl.geometryId().sensitive() / 10;
                        
                        // If the layer id of the true source link is less than
                        // the layer id of the measurement source link
                        // it means that the pivot is mismatched
                        if (trueLayerId < layerId) {
                            m_mismatch = 3;
                            norm--;
                            continue;
                        }
                        // If the layer id of the true source link is greater than
                        // the layer id of the measurement source link
                        // it means that the last measurement is missed
                        else if (trueLayerId > layerId) {
                            m_mismatch = 4;
                            trueSourceLinks.pop_back();
                            norm--;
                            continue;
                        }
                        // Do nothing if the mismatch happened in the same layer
                        // as it is a regular CKF mismatch
                        trueSourceLinks.pop_back();
                    }
                    else {
                        m_mismatch = 5;
                    }

                    // Get the true hit
                    auto trueHit = trueSsl.parameters;

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

                    // Get the distances between the true and the predicted hits
                    auto truePredictedDistance = truePredictedResidual.norm();
                    auto trueFilteredDistance = trueFilteredResidual.norm();
                    auto trueSmoothedDistance = trueSmoothedResidual.norm();

                    // Get the distances between the measurements and the predicted hits
                    auto predictedDistance = predictedResidual.norm();
                    auto filteredDistance = filteredResidual.norm();
                    auto smoothedDistance = smoothedResidual.norm();

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

                    // // KF pulls with respect to the true hits
                    // auto truePredictedPull = predictedCovariance.inverse().cwiseSqrt() * truePredictedResidual;
                    // auto trueFilteredPull = filteredCovariance.inverse().cwiseSqrt() * trueFilteredResidual;
                    // auto trueSmoothedPull = smoothedCovariance.inverse().cwiseSqrt() * trueSmoothedResidual;

                    // // KF pulls with respect to the measurements
                    // auto predictedPull = predictedCovariance.inverse().cwiseSqrt() * predictedResidual;
                    // auto filteredPull = filteredCovariance.inverse().cwiseSqrt() * filteredResidual;
                    // auto smoothedPull = smoothedCovariance.inverse().cwiseSqrt() * smoothedResidual;

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
                        TVector3(predictedResidual.x(), predictedResidual.y(), predictedResidual.z()));
                    trueFilteredResiduals.push_back(
                        TVector3(filteredResidual.x(), filteredResidual.y(), filteredResidual.z()));
                    trueSmoothedResiduals.push_back(
                        TVector3(smoothedResidual.x(), smoothedResidual.y(), smoothedResidual.z()));

                    // Store the residuals with respect to the measurements
                    predictedResiduals.push_back(
                        TVector3(predictedResidual.x(), predictedResidual.y(), predictedResidual.z()));
                    filteredResiduals.push_back(
                        TVector3(filteredResidual.x(), filteredResidual.y(), filteredResidual.z()));
                    smoothedResiduals.push_back(
                        TVector3(smoothedResidual.x(), smoothedResidual.y(), smoothedResidual.z()));

                    // Store the distances with respect to the true hits
                    truePredictedDistances.push_back(truePredictedDistance);
                    trueFilteredDistances.push_back(trueFilteredDistance);
                    trueSmoothedDistances.push_back(trueSmoothedDistance);

                    // Store the distances with respect to the measurements
                    predictedDistances.push_back(predictedDistance);
                    filteredDistances.push_back(filteredDistance);
                    smoothedDistances.push_back(smoothedDistance);

                    // // Store the pulls with respect to the true hits
                    // truePredictedPulls.push_back(
                        // TVector3(predictedPull.x(), 0, -predictedPull.y()));
                    // trueFilteredPulls.push_back(
                        // TVector3(filteredPull.x(), 0, -filteredPull.y()));
                    // trueSmoothedPulls.push_back(
                        // TVector3(smoothedPull.x(), 0, -smoothedPull.y()));

                    // // Store the pulls with respect to the measurements
                    // predictedPulls.push_back(
                        // TVector3(predictedPull.x(), 0, -predictedPull.y()));
                    // filteredPulls.push_back(
                        // TVector3(filteredPull.x(), 0, -filteredPull.y()));
                    // smoothedPulls.push_back(
                        // TVector3(smoothedPull.x(), 0, -smoothedPull.y()));
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

                // KF distances with respect to the true hits
                m_truePredictedDistances = truePredictedDistances;
                m_trueFilteredDistances = trueFilteredDistances;
                m_trueSmoothedDistances = trueSmoothedDistances;

                // KF distances with respect to the measurements
                m_predictedDistances = predictedDistances;
                m_filteredDistances = filteredDistances;
                m_smoothedDistances = smoothedDistances;

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

                m_trackId = id;

                // Matching degree
                m_matchingDegree = matchingDegree / norm;

                // Fill the tree
                m_tree->Fill();
            }

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
                m_inputTracks{this, "InputTracks"};  

        ReadDataHandle<Seeds>
            m_inputSeeds{this, "InputSeeds"};  

        std::unique_ptr<const Acts::Logger> m_logger;

        /// mutex used to protect multi-threaded reads
        std::mutex m_write_mutex;

        /// The output file
        TFile *m_file = nullptr;

        /// The output tree
        TTree *m_tree = nullptr;

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

        /// KF distances with respect to the true hits
        std::vector<double> m_truePredictedDistances;
        std::vector<double> m_trueFilteredDistances;
        std::vector<double> m_trueSmoothedDistances;

        /// KF distances with respect to the measurements
        std::vector<double> m_predictedDistances;
        std::vector<double> m_filteredDistances;
        std::vector<double> m_smoothedDistances;

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

        /// Flag indicating the number of states mismatch
        int m_mismatch;        

        /// TrackId
        int m_trackId;

        /// KF predicted momentum at the IP
        TLorentzVector m_ipMomentum;
        TVector3 m_ipMomentumError;
        TVector3 m_vertex;
        TVector3 m_vertexError;

        /// True momentum at the IP
        TLorentzVector m_ipMomentumTruth;
        TVector3 m_vertexTruth;

        /// Mutex to protect the tree filling
        std::mutex m_mutex;
};
