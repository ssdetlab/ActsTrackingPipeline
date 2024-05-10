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
        struct ROOTMeasurement {
            unsigned int id;
            float x1;
            float y1;
            float z1;
            float x4;
            float y4;
            float z4;
            float E;
        };

        /// @brief The nested configuration struct
        struct Config {
            /// The input collection
            std::string inputCollection = "Measurements";
            /// Geometry options
            const LUXEGeometry::GeometryOptions& gOpt;
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

                m_tree->Branch("id", &m_s.id);
                m_tree->Branch("x1", &m_s.x1);
                m_tree->Branch("y1", &m_s.y1);
                m_tree->Branch("z1", &m_s.z1);
                m_tree->Branch("x4", &m_s.x4);
                m_tree->Branch("y4", &m_s.y4);
                m_tree->Branch("z4", &m_s.z4);
                m_tree->Branch("E", &m_s.E);

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
                ACTS_INFO("Empty input encountered, skipping");
                return ProcessCode::SUCCESS;
            }

            auto firstHit = input.front();
            Acts::Vector2 locFirstHit(
                firstHit.truthParameters[Acts::eBoundLoc0],
                firstHit.truthParameters[Acts::eBoundLoc1]);

            auto globFirstHit = 
                m_cfg.surfaceAccessor(firstHit.sourceLink)->localToGlobal(
                    ctx.geoContext, 
                    locFirstHit, 
                    Acts::Vector3(0, 1, 0));

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

                m_s.id = firstHit.trackId;
                m_s.x1 = globFirstHit.x();
                m_s.y1 = globFirstHit.y();
                m_s.z1 = globFirstHit.z();
                m_s.x4 = globLastHit.x();
                m_s.y4 = globLastHit.y();
                m_s.z4 = globLastHit.z();
                m_s.E = std::hypot(
                    1/firstHit.truthParameters[Acts::eBoundQOverP], me);
                m_tree->Fill();
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

        /// Private access to the ROOTMeasurement struct
        ROOTMeasurement m_s;

        std::unique_ptr<const Acts::Logger> m_logger;
};
