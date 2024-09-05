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
            /// Name of the fitted track collection
            std::string inputTrackCollection;
            /// Name of the seed collection
            std::string inputSeedCollection;
            /// Name of the input tree
            std::string treeName;
            /// The names of the input files
            std::string filePath;
            /// Min number of hits
            std::uint32_t minHits = 4;
            /// Max number of hits
            std::uint32_t maxHits = 8;
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

                // Set the branches
                m_tree->Branch("trackHits", &m_trackHits, buf_size, split_lvl);

                m_tree->Branch("predictedTrackHits", &m_predictedTrackHits, buf_size, split_lvl);
                m_tree->Branch("filteredTrackHits", &m_filteredTrackHits, buf_size, split_lvl);
                m_tree->Branch("smoothedTrackHits", &m_smoothedTrackHits, buf_size, split_lvl);

                m_tree->Branch("predictedResiduals", &m_predictedResiduals, buf_size, split_lvl);
                m_tree->Branch("filteredResiduals", &m_filteredResiduals, buf_size, split_lvl);
                m_tree->Branch("smoothedResiduals", &m_smoothedResiduals, buf_size, split_lvl);

                m_tree->Branch("predictedDistances", &m_predictedDistances, buf_size, split_lvl);
                m_tree->Branch("filteredDistances", &m_filteredDistances, buf_size, split_lvl);
                m_tree->Branch("smoothedDistances", &m_smoothedDistances, buf_size, split_lvl);

                m_tree->Branch("predictedPulls", &m_predictedPulls, buf_size, split_lvl);
                m_tree->Branch("filteredPulls", &m_filteredPulls, buf_size, split_lvl);
                m_tree->Branch("smoothedPulls", &m_smoothedPulls, buf_size, split_lvl);

                m_tree->Branch("ipMomentum", &m_ipMomentum);
                m_tree->Branch("vertex", &m_vertex);

                m_tree->Branch("ipMomentumTruth", &m_ipMomentumTruth);
                m_tree->Branch("vertexTruth", &m_vertexTruth);

                m_tree->Branch("chi2", &m_chi2, "chi2/D");

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

            auto tracks = inputTracks.tracks;

            for (int idx = 0; idx < inputTracks.size(); idx++) {
                auto [id,track] = inputTracks.getByIndex(idx);

                std::vector<TVector3> trackHits;

                std::vector<TVector3> predictedTrackHits;
                std::vector<TVector3> filteredTrackHits;
                std::vector<TVector3> smoothedTrackHits;

                std::vector<TVector3> predictedResiduals;
                std::vector<TVector3> filteredResiduals;
                std::vector<TVector3> smoothedResiduals;

                std::vector<double> predictedDistances;
                std::vector<double> filteredDistances;
                std::vector<double> smoothedDistances;

                std::vector<TVector3> predictedPulls;
                std::vector<TVector3> filteredPulls;
                std::vector<TVector3> smoothedPulls;
                for (auto state : track.trackStatesReversed()) {
                    if (!state.hasProjector()) {
                        continue;
                    }
                    auto hit = state.effectiveCalibrated();

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

                    auto predictedResidual = hitGlobal - predictedHitGlobal;
                    auto filteredResidual = hitGlobal - filteredHitGlobal;
                    auto smoothedResidual = hitGlobal - smoothedHitGlobal;

                    auto predictedDistance = predictedResidual.norm();
                    auto filteredDistance = filteredResidual.norm();
                    auto smoothedDistance = smoothedResidual.norm();

                    predictedDistances.push_back(predictedDistance);
                    filteredDistances.push_back(filteredDistance);
                    smoothedDistances.push_back(smoothedDistance);

                    trackHits.push_back(
                        TVector3(hitGlobal.x(), hitGlobal.y(), hitGlobal.z()));

                    predictedTrackHits.push_back(
                        TVector3(predictedHitGlobal.x(), predictedHitGlobal.y(), predictedHitGlobal.z()));
                    filteredTrackHits.push_back(
                        TVector3(filteredHitGlobal.x(), filteredHitGlobal.y(), filteredHitGlobal.z()));
                    smoothedTrackHits.push_back(
                        TVector3(smoothedHitGlobal.x(), smoothedHitGlobal.y(), smoothedHitGlobal.z()));

                    predictedResiduals.push_back(
                        TVector3(predictedResidual.x(), predictedResidual.y(), predictedResidual.z()));
                    filteredResiduals.push_back(
                        TVector3(filteredResidual.x(), filteredResidual.y(), filteredResidual.z()));
                    smoothedResiduals.push_back(
                        TVector3(smoothedResidual.x(), smoothedResidual.y(), smoothedResidual.z()));

                    // auto predictedCovariance =
                        // state.effectiveProjector() * 
                        // state.predictedCovariance() * 
                        // state.effectiveProjector().transpose();
                    // auto predictedPull = predictedCovariance.inverse().cwiseSqrt() * predictedResidual;

                    // auto filteredCovariance =
                        // state.effectiveProjector() * 
                        // state.filteredCovariance() * 
                        // state.effectiveProjector().transpose();
                    // auto filteredPull = filteredCovariance.inverse().cwiseSqrt() * filteredResidual;

                    // auto smoothedCovariance =
                        // state.effectiveProjector() * 
                        // state.smoothedCovariance() * 
                        // state.effectiveProjector().transpose();
                    // auto smoothedPull = smoothedCovariance.inverse().cwiseSqrt() * smoothedResidual;

                    // predictedPulls.push_back(
                        // TVector3(predictedPull.x(), 0, -predictedPull.y()));
                    // filteredPulls.push_back(
                        // TVector3(filteredPull.x(), 0, -filteredPull.y()));
                    // smoothedPulls.push_back(
                        // TVector3(smoothedPull.x(), 0, -smoothedPull.y()));
                }
                m_trackHits = trackHits;
                
                m_predictedTrackHits = predictedTrackHits;
                m_filteredTrackHits = filteredTrackHits;
                m_smoothedTrackHits = smoothedTrackHits;
                
                m_predictedResiduals = predictedResiduals;
                m_filteredResiduals = filteredResiduals;
                m_smoothedResiduals = smoothedResiduals;

                m_predictedDistances = predictedDistances;
                m_filteredDistances = filteredDistances;
                m_smoothedDistances = smoothedDistances;

                m_predictedPulls = predictedPulls;
                m_filteredPulls = filteredPulls;
                m_smoothedPulls = smoothedPulls;

                m_chi2 = track.chi2();

                double me = 0.511 * Acts::UnitConstants::MeV;
                Acts::Vector3 pVec = track.momentum();
                double pMag = pVec.norm();
                m_ipMomentum.SetPxPyPzE(pVec.x(), pVec.y(), pVec.z(), std::hypot(pMag, me));
                
                Acts::Vector3 vertex = {track.loc0(), 0, -track.loc1()};
                m_vertex = TVector3(vertex.x(), vertex.y(), vertex.z());

                for (auto& seed : inputSeeds) {
                    if (seed.trackId == id) {
                        Acts::Vector3 pVecTruth = seed.ipParameters.momentum();
                        double pMagTruth = pVecTruth.norm();
                        m_ipMomentumTruth.SetPxPyPzE(
                            pVecTruth.x(), pVecTruth.y(), pVecTruth.z(), std::hypot(pMagTruth, me));
                        
                        Acts::Vector3 vertexTruth = {
                            seed.ipParameters.position().x(),
                            seed.ipParameters.position().y(),
                            seed.ipParameters.position().z()};
                        m_vertexTruth = TVector3(vertexTruth.x(), vertexTruth.y(), vertexTruth.z());

                        break;
                    }
                }

                {
                    std::lock_guard<std::mutex> lock(m_mutex);
                    m_tree->Fill();
                }
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
        /// The exausitive list of columns
        std::vector<TVector3> m_trackHits;
        
        std::vector<TVector3> m_predictedTrackHits;
        std::vector<TVector3> m_filteredTrackHits;
        std::vector<TVector3> m_smoothedTrackHits;

        std::vector<TVector3> m_predictedResiduals;
        std::vector<TVector3> m_filteredResiduals;
        std::vector<TVector3> m_smoothedResiduals;

        std::vector<double> m_predictedDistances;
        std::vector<double> m_filteredDistances;
        std::vector<double> m_smoothedDistances;

        std::vector<TVector3> m_predictedPulls;
        std::vector<TVector3> m_filteredPulls;
        std::vector<TVector3> m_smoothedPulls;

        double m_chi2;

        TLorentzVector m_ipMomentum;
        TVector3 m_vertex;

        TLorentzVector m_ipMomentumTruth;
        TVector3 m_vertexTruth;

        std::mutex m_mutex;

};
