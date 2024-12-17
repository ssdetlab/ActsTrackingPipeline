#include "TrackingPipeline/Io/RootTrackLookupValidationWriter.hpp"

RootTrackLookupValidationWriter::RootTrackLookupValidationWriter(
    const Config& config, Acts::Logging::Level level)
        : m_cfg(config), 
        m_logger(Acts::getDefaultLogger(name(), level)) {
            m_outputFile = 
                TFile::Open(m_cfg.path.c_str(), "RECREATE");
            m_outputTree = 
                new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());

            m_outputTree->Branch("idx", &m_event.idx);
        
            m_outputTree->Branch("true_vertexX", &m_event.true_vertexX);
            m_outputTree->Branch("true_vertexY", &m_event.true_vertexY);
            m_outputTree->Branch("true_vertexZ", &m_event.true_vertexZ);
        
            m_outputTree->Branch("est_vertexX", &m_event.est_vertexX);
            m_outputTree->Branch("est_vertexY", &m_event.est_vertexY);
            m_outputTree->Branch("est_vertexZ", &m_event.est_vertexZ);
        
            m_outputTree->Branch("true_refX", &m_event.true_refX);
            m_outputTree->Branch("true_refY", &m_event.true_refY);
            m_outputTree->Branch("true_refZ", &m_event.true_refZ);
        
            m_outputTree->Branch("est_refX", &m_event.est_refX);
            m_outputTree->Branch("est_refY", &m_event.est_refY);
            m_outputTree->Branch("est_refZ", &m_event.est_refZ);
        
            m_outputTree->Branch("true_IpPx", &m_event.true_IpPx);
            m_outputTree->Branch("true_IpPy", &m_event.true_IpPy);
            m_outputTree->Branch("true_IpPz", &m_event.true_IpPz);
            m_outputTree->Branch("true_IpE", &m_event.true_IpE);
        
            m_outputTree->Branch("est_IpPx", &m_event.est_IpPx);
            m_outputTree->Branch("est_IpPy", &m_event.est_IpPy);
            m_outputTree->Branch("est_IpPz", &m_event.est_IpPz);
            m_outputTree->Branch("est_IpE", &m_event.est_IpE);
        
            m_outputTree->Branch("true_refPx", &m_event.true_refPx);
            m_outputTree->Branch("true_refPy", &m_event.true_refPy);
            m_outputTree->Branch("true_refPz", &m_event.true_refPz);
            m_outputTree->Branch("true_refE", &m_event.true_refE);
        
            m_outputTree->Branch("est_refPx", &m_event.est_refPx);
            m_outputTree->Branch("est_refPy", &m_event.est_refPy);
            m_outputTree->Branch("est_refPz", &m_event.est_refPz);
            m_outputTree->Branch("est_refE", &m_event.est_refE);
        
            m_inputIpPars.initialize(m_cfg.inputIpPars);
            m_inputRefLayerPars.initialize(m_cfg.inputRefLayerPars);
            m_inputIpParsEst.initialize(m_cfg.inputIpParsEst);
            m_inputRefLayerParsEst.initialize(m_cfg.inputRefLayerParsEst);
}

ProcessCode RootTrackLookupValidationWriter::finalize() {
    m_outputFile->Write();
    m_outputFile->Close();

    return ProcessCode::SUCCESS;
}

ProcessCode 
RootTrackLookupValidationWriter::write(const AlgorithmContext& ctx) {
    auto ipPars = m_inputIpPars(ctx);
    auto refLayerPars = m_inputRefLayerPars(ctx);
    auto ipParsEst = m_inputIpParsEst(ctx);
    auto refLayerParsEst = m_inputRefLayerParsEst(ctx);

    std::lock_guard<std::mutex> lock(m_writeMutex);

    for (int i = 0; i < ipPars.size(); i++) {
        m_event.idx = i;

        m_event.true_vertexX = ipPars[i].fourPosition().x();
        m_event.true_vertexY = ipPars[i].fourPosition().y();
        m_event.true_vertexZ = ipPars[i].fourPosition().z();
    
        m_event.est_vertexX = ipParsEst[i].fourPosition().x();
        m_event.est_vertexY = ipParsEst[i].fourPosition().y();
        m_event.est_vertexZ = ipParsEst[i].fourPosition().z();
    
        m_event.true_refX = refLayerPars[i].fourPosition().x();
        m_event.true_refY = refLayerPars[i].fourPosition().y();
        m_event.true_refZ = refLayerPars[i].fourPosition().z();
    
        m_event.est_refX = refLayerParsEst[i].fourPosition().x();
        m_event.est_refY = refLayerParsEst[i].fourPosition().y();
        m_event.est_refZ = refLayerParsEst[i].fourPosition().z();
    
        m_event.true_IpPx = ipPars[i].momentum().x();
        m_event.true_IpPy = ipPars[i].momentum().y();
        m_event.true_IpPz = ipPars[i].momentum().z();
        m_event.true_IpE = ipPars[i].absoluteMomentum();
    
        m_event.est_IpPx = ipParsEst[i].momentum().x();
        m_event.est_IpPy = ipParsEst[i].momentum().y();
        m_event.est_IpPz = ipParsEst[i].momentum().z();
        m_event.est_IpE = ipParsEst[i].absoluteMomentum();
    
        m_event.true_refPx = refLayerPars[i].momentum().x();
        m_event.true_refPy = refLayerPars[i].momentum().y();
        m_event.true_refPz = refLayerPars[i].momentum().z();
        m_event.true_refE = refLayerPars[i].absoluteMomentum();
    
        m_event.est_refPx = refLayerParsEst[i].momentum().x();
        m_event.est_refPy = refLayerParsEst[i].momentum().y();
        m_event.est_refPz = refLayerParsEst[i].momentum().z();
        m_event.est_refE = refLayerParsEst[i].absoluteMomentum();
    
        m_outputTree->Fill();
    }

    // Write the track parameters to the output file
    return ProcessCode::SUCCESS;
}
