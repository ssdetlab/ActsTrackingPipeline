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
            /// name of the fitted track collection
            std::string inputTrackCollection;
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

                m_tree->Branch("predictedPulls", &m_predictedPulls, buf_size, split_lvl);
                m_tree->Branch("filteredPulls", &m_filteredPulls, buf_size, split_lvl);
                m_tree->Branch("smoothedPulls", &m_smoothedPulls, buf_size, split_lvl);

                m_inputTracks.initialize(m_cfg.inputTrackCollection);
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

            auto tracks = inputTracks.tracks;
            auto ids = inputTracks.trackIds;

            for (int idx = 0; idx < inputTracks.size(); idx++) {
                auto [id,track] = inputTracks.getByIndex(idx);
                std::cout << "WRITE ID = " << id << std::endl;

                std::vector<TVector3> trackHits;

                std::vector<TVector3> predictedTrackHits;
                std::vector<TVector3> filteredTrackHits;
                std::vector<TVector3> smoothedTrackHits;

                std::vector<TVector3> predictedResiduals;
                std::vector<TVector3> filteredResiduals;
                std::vector<TVector3> smoothedResiduals;

                std::vector<TVector3> predictedPulls;
                std::vector<TVector3> filteredPulls;
                std::vector<TVector3> smoothedPulls;
                for (auto state : track.trackStatesReversed()) {
                    if (!state.hasReferenceSurface() || !state.hasProjector()) {
                        continue;
                    }
                    auto hit = state.effectiveCalibrated();

                    auto predictedHit = state.effectiveProjector() * state.predicted();
                    auto filteredHit = state.effectiveProjector() * state.filtered();
                    auto smoothedHit = state.effectiveProjector() * state.smoothed();

                    auto predictedResidual = (hit - predictedHit).eval();
                    auto filteredResidual = (hit - filteredHit).eval();
                    auto smoothedResidual = (hit - smoothedHit).eval();

                    auto predictedDistance = predictedResidual.norm();
                    auto filteredDistance = filteredResidual.norm();
                    auto smoothedDistance = smoothedResidual.norm();

                    trackHits.push_back(TVector3(hit[0], hit[1], 0));

                    predictedTrackHits.push_back(TVector3(predictedHit[0], predictedHit[1], 0));
                    filteredTrackHits.push_back(TVector3(filteredHit[0], filteredHit[1], 0));
                    smoothedTrackHits.push_back(TVector3(smoothedHit[0], smoothedHit[1], 0));

                    predictedResiduals.push_back(TVector3(predictedResidual[0], predictedResidual[1], 0));
                    filteredResiduals.push_back(TVector3(filteredResidual[0], filteredResidual[1], 0));
                    smoothedResiduals.push_back(TVector3(smoothedResidual[0], smoothedResidual[1], 0));

                    // auto predictedCovariance =
                        // state.effectiveProjector() * 
                        // state.predictedCovariance() * 
                        // state.effectiveProjector().transpose();
                    // auto predictedPull = predictedResidual.transpose() * predictedCovariance.inverse() * predictedResidual;

                    // auto filteredCovariance =
                        // state.effectiveProjector() * 
                        // state.filteredCovariance() * 
                        // state.effectiveProjector().transpose();
                    // auto filteredPull = filteredResidual.transpose() * filteredCovariance.inverse() * filteredResidual;

                    // auto smoothedCovariance =
                        // state.effectiveProjector() * 
                        // state.smoothedCovariance() * 
                        // state.effectiveProjector().transpose();
                    // auto smoothedPull = smoothedResidual.transpose() * smoothedCovariance.inverse() * smoothedResidual;

                    // predictedPulls.push_back(TVector3(predictedPull[0], predictedPull[1], 0));
                    // filteredPulls.push_back(TVector3(filteredPull[0], filteredPull[1], 0));
                    // smoothedPulls.push_back(TVector3(smoothedPull[0], smoothedPull[1], 0));
                }
                m_trackHits = trackHits;
                m_predictedTrackHits = predictedTrackHits;
                m_filteredTrackHits = filteredTrackHits;
                m_smoothedTrackHits = smoothedTrackHits;
                m_predictedResiduals = predictedResiduals;
                m_filteredResiduals = filteredResiduals;
                m_smoothedResiduals = smoothedResiduals;
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

        std::vector<TVector3> m_predictedPulls;
        std::vector<TVector3> m_filteredPulls;
        std::vector<TVector3> m_smoothedPulls;


        // std::vector<TMatrix> m_predictedCovariances;
        // std::vector<TMatrix> m_filteredCovariances;
        // std::vector<TMatrix> m_smoothedCovariances;
};
