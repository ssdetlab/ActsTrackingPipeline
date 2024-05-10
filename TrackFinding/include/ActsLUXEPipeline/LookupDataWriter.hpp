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

class LookupDataWriter : public IWriter {
    public:
        struct ROOTMeasurement {
            unsigned int id;
            float x1;
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
            /// Detector
            const std::shared_ptr<const Acts::Experimental::Detector>& detector;
            /// The output file name
            std::string filePath = "seed-lookup.root";
            /// The output tree name
            std::string treeName = "seed-lookup";
        };


        /// @brief Constructor
        LookupDataWriter(const Config& config, Acts::Logging::Level level)
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
                m_tree->Branch("z1", &m_s.z1);
                m_tree->Branch("x4", &m_s.x4);
                m_tree->Branch("y4", &m_s.y4);
                m_tree->Branch("z4", &m_s.z4);
                m_tree->Branch("E", &m_s.E);

                m_inputMeasurements.initialize(m_cfg.inputCollection);
        }

        LookupDataWriter(const LookupDataWriter &) = delete;
        LookupDataWriter(const LookupDataWriter &&) = delete;

        ~LookupDataWriter() {
            if (m_file) {
                m_file->Write();
                m_file->Close();
                delete m_file;
            }
        }

        /// Writer name() method
        std::string name() const { return "LookupDataWriter"; }

        /// @brief The execute method        
        ProcessCode write(const AlgorithmContext &ctx) override {
            // Get the input measurements
            // from the context
            auto input = m_inputMeasurements(ctx);

            SimpleSourceLink::SurfaceAccessor SA{*m_cfg.detector};

            if (input.empty()) {
                ACTS_INFO("Empty input encountered, skipping");
                return ProcessCode::SUCCESS;
            }

            Acts::Vector2 params(
                input[0].truthParameters[0],
                input[0].truthParameters[1]);

            auto globalPosition = 
                SA(input[0].sourceLink)->localToGlobal(
                    ctx.geoContext, 
                    params, 
                    Acts::Vector3(0, 1, 0));

            auto me = 0.511 * Acts::UnitConstants::MeV;

            auto size = input.size();
            if ((globalPosition[1] < m_cfg.gOpt.staveZ.at(3) && 
                globalPosition[0] <= m_cfg.gOpt.chipXOdd.at(0)) ||
                    globalPosition[1] < m_cfg.gOpt.staveZ.at(0)) {
                        Acts::Vector2 paramsL(
                            input[size-1].truthParameters[0],
                            input[size-1].truthParameters[1]);

                        auto globalPositionL = 
                            SA(input[size-1].sourceLink)->localToGlobal(
                                ctx.geoContext, 
                                paramsL, 
                                Acts::Vector3(0, 1, 0));

                        m_s.id = input[0].trackId;
                        m_s.x1 = globalPosition[0];
                        m_s.z1 = globalPosition[2];
                        m_s.x4 = globalPositionL[0];
                        m_s.y4 = globalPositionL[1];
                        m_s.z4 = globalPositionL[2];
                        m_s.E = std::hypot(
                            std::pow(input[0].truthParameters[4],-1), me);
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
