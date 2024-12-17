#pragma once

#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h" 

using namespace Acts::UnitLiterals;

using TrackID = std::tuple<std::int32_t, std::int32_t, std::int32_t>;

/// @brief Writer to store fitted track data in
/// ROOT file
///
/// Writer that accepts fitted track data from KF 
/// derives the basic performance metrics, such as
/// chi2 and residuals, and stores them in a ROOT file.
///
/// @note Assumes that the tracks are simulated and 
/// the truth information is available
class RootFittedSimTrackWriter : public IWriter {
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
            /// Target size of the true track
            std::size_t targetTrueTrackSize;
        };

        RootFittedSimTrackWriter(const RootFittedSimTrackWriter &) = delete;
        RootFittedSimTrackWriter(const RootFittedSimTrackWriter &&) = delete;
    
        /// @brief Constructor
        ///
        /// @param config The Configuration struct
        RootFittedSimTrackWriter(const Config &config, Acts::Logging::Level level);

        /// @brief Finalize method
        ProcessCode finalize() override; 
    
        /// Writer name() method
        std::string name() const override { return "RootFittedTrackWriter"; }
    
        /// Write out data to the input stream
        ProcessCode write(const AlgorithmContext &ctx) override; 

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

        /// KF pulls with respect to the true hits
        std::vector<TVector3> m_truePredictedPulls;
        std::vector<TVector3> m_trueFilteredPulls;
        std::vector<TVector3> m_trueSmoothedPulls;

        /// KF pulls with respect to the measurements
        std::vector<TVector3> m_predictedPulls;
        std::vector<TVector3> m_filteredPulls;
        std::vector<TVector3> m_smoothedPulls;

        /// Chi2 of the track
        /// with respect to the 
        /// true hit position
        double m_trueChi2;

        /// Chi2 of the track
        /// with respect ot the 
        /// measurement
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

        /// Mutex to protect the tree filling
        std::mutex m_mutex;
};
