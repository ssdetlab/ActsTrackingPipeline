#pragma once

#include "Acts/EventData/SourceLink.hpp"

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"
#include "ActsLUXEPipeline/LUXEDataContainers.hpp"
#include "ActsLUXEPipeline/LUXEGeometryConstraints.hpp"
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

/// @brief The ideal seeder for the LUXE simulation
/// takes the the SimMeasurements and converts them 
/// into seeds

class LookupTableMaker : public IAlgorithm {
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
            /// filename
            std::string filename = "output_file.root";

            TFile* file;
        };


        /// @brief Constructor
        LookupTableMaker(Config config, Acts::Logging::Level level)
            : IAlgorithm("LookupTableMaker", level),
            m_cfg(std::move(config)) {
                m_inputMeasurements.initialize(m_cfg.inputCollection);
                ACTS_VERBOSE("Constructing lookup table...");
                m_cfg.file = new TFile(&m_cfg.filename[0],"recreate");
        }

        ~LookupTableMaker() = default;

        /// @brief The execute method        
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            using namespace Acts::UnitLiterals;

            // Get the input measurements
            // from the context
            auto input = m_inputMeasurements(ctx);


            TTree* tree = dynamic_cast<TTree*>(m_cfg.file->Get("basic_tree"));
            ROOTMeasurement s;
            if (!tree) {
                tree = new TTree("basic_tree", "TruthParameters");
                tree->Branch("id", &s.id);
                tree->Branch("x1", &s.x1);
                tree->Branch("z1", &s.z1);
                tree->Branch("x4", &s.x4);
                tree->Branch("y4", &s.y4);
                tree->Branch("z4", &s.z4);
                tree->Branch("E", &s.E);
            }

            SimpleSourceLink::SurfaceAccessor SA{*m_cfg.detector};

            if (input.size()==0) {
                ACTS_INFO("Empty input encountered, skipping");
                return ProcessCode::SUCCESS;
            }

            Acts::Vector2 params{input[0].truthParameters[0],input[0].truthParameters[1]};
            auto globalPosition = SA(input[0].sourceLink)->
                    localToGlobal(ctx.geoContext, params, Acts::Vector3{0,1,0});
            auto size = input.size();
            if ((globalPosition[1]<m_cfg.gOpt.staveZ.at(3) && globalPosition[0]<=m_cfg.gOpt.chipXOdd.at(0)) ||
                    globalPosition[1]<m_cfg.gOpt.staveZ.at(0)) {
                Acts::Vector2 paramsL{input[size-1].truthParameters[0],input[size-1].truthParameters[1]};
                auto globalPositionL = SA(input[size-1].sourceLink)->
                        localToGlobal(ctx.geoContext, paramsL, Acts::Vector3{0,1,0});;
                s.id = input[0].trackId;
                s.x1 = globalPosition[0];
                s.z1 = globalPosition[2];
                s.E = std::hypot(std::pow(input[0].truthParameters[4],-1),0.000511);
                s.x4 = globalPositionL[0];
                s.y4 = globalPositionL[1];
                s.z4 = globalPositionL[2];
                tree->Fill();
            }
            m_cfg.file->Write(0,TObject::kOverwrite);
            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;

        ReadDataHandle<LUXEDataContainer::SimMeasurements> m_inputMeasurements
            {this, "InputMeasurements"};
};
