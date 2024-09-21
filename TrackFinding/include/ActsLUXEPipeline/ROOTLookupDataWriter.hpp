#pragma once

#include "ActsLUXEPipeline/IWriter.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"
#include "ActsLUXEPipeline/LUXEGeometryConstraints.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Definitions/Units.hpp" 

#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TBrowser.h>
#include <TStyle.h>
#include <TLatex.h>

#include <vector>
#include <iostream>

class ROOTLookupDataWriter : public IWriter {
    public:
        /// @brief The ROOT lookup data struct
        /// for seeding Lookup table creation
        struct ROOTLookupData {
            unsigned int id;
            float xFirst;
            float yFirst;
            float zFirst;
            float xLast;
            float yLast;
            float zLast;
            float ipPx;
            float ipPy;
            float ipPz;
            float E;
        };

        /// @brief The nested configuration struct
        struct Config {
            /// The input collection
            std::string inputCollection = "Measurements";
            /// Surface accessor
            Acts::SourceLinkSurfaceAccessor surfaceAccessor;
            /// First layer extent
            Acts::Extent firstLayerExtent;
            /// Last layer extent
            Acts::Extent lastLayerExtent;
            /// The output file name
            std::string filePath = "seed-lookup.root";
            /// The output tree name
            std::string treeName = "seed-lookup";
        };


        /// @brief Constructor
        ROOTLookupDataWriter(const Config& config, Acts::Logging::Level level)
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

                m_tree->Branch("id", &m_ld.id);
                m_tree->Branch("xFirst", &m_ld.xFirst);
                m_tree->Branch("yFirst", &m_ld.yFirst);
                m_tree->Branch("zFirst", &m_ld.zFirst);
                m_tree->Branch("xLast", &m_ld.xLast);
                m_tree->Branch("yLast", &m_ld.yLast);
                m_tree->Branch("zLast", &m_ld.zLast);
                m_tree->Branch("ipPx", &m_ld.ipPx);
                m_tree->Branch("ipPy", &m_ld.ipPy);
                m_tree->Branch("ipPz", &m_ld.ipPz);
                m_tree->Branch("E", &m_ld.E);

                m_inputMeasurements.initialize(m_cfg.inputCollection);
        }

        ROOTLookupDataWriter(const ROOTLookupDataWriter &) = delete;
        ROOTLookupDataWriter(const ROOTLookupDataWriter &&) = delete;

        ~ROOTLookupDataWriter() {
            if (m_file) {
                m_file->Write();
                m_file->Close();
                delete m_file;
            }
        }

        /// Writer name() method
        std::string name() const { return "ROOTLookupDataWriter"; }

        /// @brief The execute method        
        ProcessCode write(const AlgorithmContext &ctx) override {
            // Get the input measurements
            // from the context
            auto input = m_inputMeasurements(ctx);

            if (input.empty()) {
                return ProcessCode::SUCCESS;
            }

            // Sort the input by y coordinate
            std::sort(input.begin(), input.end(), 
                [&](const SimMeasurement& a, const SimMeasurement& b) {
                    Acts::Vector2 aHitLoc(
                        a.truthParameters[Acts::eBoundLoc0],
                        a.truthParameters[Acts::eBoundLoc1]);

                    Acts::Vector2 bHitLoc(
                        b.truthParameters[Acts::eBoundLoc0],
                        b.truthParameters[Acts::eBoundLoc1]);

                    auto aHitGlob = 
                        m_cfg.surfaceAccessor(a.sourceLink)->localToGlobal(
                            ctx.geoContext, 
                            aHitLoc, 
                            Acts::Vector3(0, 1, 0));
                    
                    auto bHitGlob =
                        m_cfg.surfaceAccessor(b.sourceLink)->localToGlobal(
                            ctx.geoContext, 
                            bHitLoc, 
                            Acts::Vector3(0, 1, 0));

                    return aHitGlob.y() < bHitGlob.y();
                });

            auto firstHit = input.front();
            Acts::Vector2 locFirstHit(
                firstHit.truthParameters[Acts::eBoundLoc0],
                firstHit.truthParameters[Acts::eBoundLoc1]);

            auto globFirstHit = 
                m_cfg.surfaceAccessor(firstHit.sourceLink)->localToGlobal(
                    ctx.geoContext, 
                    locFirstHit, 
                    Acts::Vector3(0, 1, 0));

            // Store the lookup data
            auto me = 0.511 * Acts::UnitConstants::MeV;
            if (m_cfg.firstLayerExtent.contains(globFirstHit)) {
                auto lastHit = input.back();

                Acts::Vector2 locLastHit(
                    lastHit.truthParameters[Acts::eBoundLoc0],
                    lastHit.truthParameters[Acts::eBoundLoc1]);

                auto globLastHit = 
                    m_cfg.surfaceAccessor(lastHit.sourceLink)->localToGlobal(
                        ctx.geoContext, 
                        locLastHit, 
                        Acts::Vector3(0, 1, 0));

                {
                    std::lock_guard<std::mutex> lock(m_mutex);

                    m_ld.id = firstHit.trackId;
                    m_ld.xFirst = globFirstHit.x();
                    m_ld.yFirst = globFirstHit.y();
                    m_ld.zFirst = globFirstHit.z();
                    m_ld.xLast = globLastHit.x();
                    m_ld.yLast = globLastHit.y();
                    m_ld.zLast = globLastHit.z();
                    m_ld.ipPx = firstHit.ipParameters.momentum().x();
                    m_ld.ipPy = firstHit.ipParameters.momentum().y();
                    m_ld.ipPz = firstHit.ipParameters.momentum().z();
                    m_ld.E = std::hypot(
                        1/firstHit.truthParameters[Acts::eBoundQOverP], me);
                    m_tree->Fill();
                }
            }

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        /// The config class
        Config m_cfg;

        /// Private access to the logging instance
        const Acts::Logger &logger() const { return *m_logger; }

        ReadDataHandle<SimMeasurements> m_inputMeasurements
            {this, "InputMeasurements"};

        /// The output file
        TFile *m_file = nullptr;

        /// The output tree
        TTree *m_tree = nullptr;

        /// Private access to the ROOTLookupData struct
        ROOTLookupData m_ld;

        std::unique_ptr<const Acts::Logger> m_logger;

        std::mutex m_mutex;
};
