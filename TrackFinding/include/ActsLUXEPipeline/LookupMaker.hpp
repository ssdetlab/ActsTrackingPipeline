#pragma once

#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <TH2D.h>
#include <TInterpreter.h>
#include <TROOT.h>

#include <vector>
#include <iostream>

class LookupMaker {
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
            std::string filename = "output_file.root";
            /// Geometry options
        };


        /// @brief Constructor
        LookupMaker(Config config)
            : m_cfg(std::move(config)) {
//                m_cfg.file = new TFile(&m_cfg.filename[0],"recreate");
        }

        ~LookupMaker() = default;

        std::unordered_map<float, float> makeLookUp(TH2D* histogram, std::string filename) const {
            std::unordered_map<float, float> lookupTable;
            // Loop over bins to build the lookup table
            for (int binx = 1; binx <= histogram->GetNbinsX(); ++binx) {
                float x = histogram->GetXaxis()->GetBinCenter(binx);
                float cumSum = 0.0;
                std::vector<float> values;
                for (int biny = 1; biny <= histogram->GetNbinsY(); ++biny) {
                    cumSum += histogram->GetBinContent(binx, biny);
                    values.push_back(cumSum);
                }
                float totalSum = cumSum;
                auto it = std::lower_bound(values.begin(), values.end(), cumSum*0.5);

                float median = histogram->GetYaxis()->GetBinCenter(std::distance(values.begin(), it));
                lookupTable[x] = median;
            }
            FILE *lookupFile = fopen(&filename[0], "w");
            if (lookupFile) {
                for (const auto &entry : lookupTable) {
                    fprintf(lookupFile, "%f %f\n", entry.first, entry.second);
                }
                fclose(lookupFile);
            } else {
                std::cerr << "Failed to open lookup table file for writing." << std::endl;
            }
            return lookupTable;
        }

        /// @brief The execute method        
        void execute() const {
            TROOT root("root", "ROOT program");
            // Get the input measurements
            // from the context

            gInterpreter->Declare("");
            TFile *file = new TFile(&m_cfg.filename[0]);

            // Retrieve the tree from the file
            TTree *tree;
            file->GetObject("basic_tree", tree);

            Float_t x1_data, z1_data, x4_data, y4_data, z4_data, E_data;
            tree->SetBranchAddress("x1", &x1_data);
            tree->SetBranchAddress("z1", &z1_data);
            tree->SetBranchAddress("x4", &x4_data);
            tree->SetBranchAddress("y4", &y4_data);
            tree->SetBranchAddress("z4", &z4_data);
            tree->SetBranchAddress("E", &E_data);

            // Create a 2D histogram
            TH2D *EX1histogram = new TH2D("EX1histogram", "", 1500, 50, 560, 400, 0.5, 16.5); // Adjust binning and range as needed

            // Fill the histogram from the tree
            Long64_t nEntries = tree->GetEntries();
            for (Long64_t i = 0; i < nEntries; ++i) {
                tree->GetEntry(i);
                EX1histogram->Fill(x1_data, E_data);
            }

            TH2D *X1X4histogram = new TH2D("X1X4histogram", "", 3000, 50, 480, 3000, 50, 560); // Adjust binning and range as needed

            // Fill histogram
            for (Long64_t i = 0; i < nEntries; ++i) {
                tree->GetEntry(i);
                X1X4histogram->Fill(x1_data, x4_data);
            }

            TH2D *X1Y4histogram = new TH2D("X1Y4histogram", "", 2500, 50, 560, 13, 4249.5125, 4262.5125); // Adjust binning and range as needed
            // Fill histogram
            for (Long64_t i = 0; i < nEntries; ++i) {
                tree->GetEntry(i);
                X1Y4histogram->Fill(x1_data, y4_data);
            }

            TH2D *Z1Z4histogram = new TH2D("Z1Z4histogram", "", 200, -6, 6, 200, -6, 6); // Adjust binning and range as needed
            // Fill histogram
            for (Long64_t i = 0; i < nEntries; ++i) {
                tree->GetEntry(i);
                Z1Z4histogram->Fill(z1_data, z4_data);
            }

            auto EX1_lookup = makeLookUp(EX1histogram,"EX1_lookup_table.txt");
            auto X1X4_lookup = makeLookUp(X1X4histogram,"X1X4_lookup_table.txt");
            auto X1Y4_lookup = makeLookUp(X1Y4histogram,"X1Y4_lookup_table.txt");
            auto Z1Z4_lookup = makeLookUp(Z1Z4histogram,"Z1Z4_lookup_table.txt");
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;
};
