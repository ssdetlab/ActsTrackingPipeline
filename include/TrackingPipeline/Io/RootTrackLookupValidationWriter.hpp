#pragma once

#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

#include "Acts/EventData/TrackParameters.hpp"

#include <mutex>

#include "TFile.h"
#include "TTree.h"

class RootTrackLookupValidationWriter final : public IWriter {
    public:
        struct EventData {
            int idx;
        
            double true_vertexX;
            double true_vertexY;
            double true_vertexZ;
        
            double est_vertexX;
            double est_vertexY;
            double est_vertexZ;
        
            double true_refX;
            double true_refY;
            double true_refZ;
        
            double est_refX;
            double est_refY;
            double est_refZ;
        
            double true_IpPx;
            double true_IpPy;
            double true_IpPz;
            double true_IpE;
        
            double est_IpPx;
            double est_IpPy;
            double est_IpPz;
            double est_IpE;
        
            double true_refPx;
            double true_refPy;
            double true_refPz;
            double true_refE;
        
            double est_refPx;
            double est_refPy;
            double est_refPz;
            double est_refE;
        };

        struct Config {
            /// Input IpPars collection
            std::string inputIpPars = "inputIpPars";
            /// Input RefLayerPars collection
            std::string inputRefLayerPars = "inputRefLayerPars";
            /// Input IpParsEst collection
            std::string inputIpParsEst = "inputIpParsEst";
            /// Input RefLayerParsEst collection
            std::string inputRefLayerParsEst = "inputRefLayerParsEst";
            /// Output file name
            std::string path = "trackParamsValidation.root";
            /// Output tree name
            std::string treeName = "track-params";
        };
        
        RootTrackLookupValidationWriter(
            const Config& config, Acts::Logging::Level level);

        ProcessCode finalize() override;

        ProcessCode write(const AlgorithmContext& ctx) override;

        /// The algorithm name.
        std::string name() const final { return "RootTrackParamsValidationWriter"; }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
        
        const Acts::Logger& logger() const { return *m_logger; }

    private:
        Config m_cfg;
        
        TFile* m_outputFile = nullptr;
        TTree* m_outputTree = nullptr;
        
        EventData m_event;
        
        std::mutex m_writeMutex;
        
        std::unique_ptr<const Acts::Logger> m_logger;
        
        ReadDataHandle<std::vector<Acts::CurvilinearTrackParameters>> m_inputIpPars{
            this, "InputIpPars"};
        ReadDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
            m_inputRefLayerPars{this, "InputRefLayerPars"};
        
        ReadDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
            m_inputIpParsEst{this, "InputIpParsEst"};
        ReadDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
            m_inputRefLayerParsEst{this, "InputRefLayerParsEst"};
};
