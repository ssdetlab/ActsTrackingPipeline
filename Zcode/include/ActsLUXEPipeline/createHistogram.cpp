#pragma once

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>

#include <TH3F.h>
#include <TCanvas.h>
#include <TBrowser.h>
#include <vector>
#include <iostream>

TGraph* CalculatePercentileCurve(TH2F *histogram, float percentile) {

    int numXBins = histogram->GetNbinsX();
    TGraph *graph = new TGraph(numXBins);

    for (int binX = 1; binX <= numXBins; ++binX) {
        std::vector<float> values;
        float cumSum = 0.0;
        for (int binY = 1; binY <= histogram->GetNbinsY(); ++binY) {
            float yValue = histogram->GetBinContent(binX, binY);
            cumSum+=yValue;
//            values.push_back(histogram->GetYaxis()->GetBinCenter(binY));
            values.push_back(cumSum);
        }

        float totalSum = cumSum;
        float percentileSum = percentile / 100.0 * totalSum;

        // Find the index where the cumulative sum equals or exceeds the percentile
        auto it = std::lower_bound(values.begin(), values.end(), percentileSum);

        // Calculate the y-value corresponding to the percentile
        float percentileValue = histogram->GetYaxis()->GetBinCenter(std::distance(values.begin(), it));

        // Set point in the graph
        graph->SetPoint(binX - 1, histogram->GetXaxis()->GetBinCenter(binX), percentileValue);
    }
    return graph;
}

std::unordered_map<float, float> makeLookUp(TH2F* histogram, std::string filename,
                                            bool yLookup = false) {

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
        if (yLookup) {
            if (median<4255) { // avg(layer6,layer7)
                median = 4250.0125; // gOpt.layerZ.at(7)
            } else {
                median = 4262.0125; // gOpt.layerZ.at(6)
            }
        }
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

//std::unordered_map<float, std::vector<float>> make3DLookUp(TH3F* histogram, std::string filename) {
//
//    std::unordered_map<float, std::vector<float>> lookupTable;
//    // Loop over bins to build the lookup table
//    for (int binx = 1; binx <= histogram->GetNbinsX(); ++binx) {
//        float x = histogram->GetXaxis()->GetBinCenter(binx);
//        std::vector<float> y_values;
//        std::vector<float> z_values(static_cast<size_t>(histogram->GetNbinsZ()));
//        float cumSum=0;
//        float ZSum;
//        for (int biny = 1; biny <= histogram->GetNbinsY(); ++biny) {
//            ZSum=0;
//            for (int binz = 1; binz <= histogram->GetNbinsZ(); ++binz) {
//                cumSum += histogram->GetBinContent(binx, biny, binz);
//                ZSum += histogram->GetBinContent(binx, biny, binz);
//                z_values[binz-1] = ZSum;
//            }
//            y_values.push_back(cumSum);
//        }
//        float totalSum = cumSum;
//        auto y_it = std::lower_bound(y_values.begin(), y_values.end(), totalSum*0.5);
//        auto z_it = std::lower_bound(z_values.begin(), z_values.end(), totalSum*0.5);
//
//        float Y_median = histogram->GetYaxis()->GetBinCenter(std::distance(y_values.begin(), y_it));
//        float Z_median = histogram->GetZaxis()->GetBinCenter(std::distance(z_values.begin(), z_it));
//        std::vector<float> v{Y_median,Z_median};
//        lookupTable[x] = v;
//    }
//    FILE *lookupFile = fopen(&filename[0], "w");
//    if (lookupFile) {
//        for (const auto &entry : lookupTable) {
//            fprintf(lookupFile, "%f %f %f\n", entry.first, entry.second[0], entry.second[1]);
//        }
//        fclose(lookupFile);
//    } else {
//        std::cerr << "Failed to open lookup table file for writing." << std::endl;
//    }
//    return lookupTable;
//}

void createHistogram() {
    TFile *file = new TFile("Zroot_files/hist_data_o.root");

    // Retrieve the tree from the file
    TTree *tree;
    file->GetObject("tree", tree);

    Float_t x1_data, z1_data, x4_data, y4_data, z4_data, E_data;
    tree->SetBranchAddress("x1", &x1_data);
    tree->SetBranchAddress("z1", &z1_data);
    tree->SetBranchAddress("x4", &x4_data);
    tree->SetBranchAddress("y4", &y4_data);
    tree->SetBranchAddress("z4", &z4_data);
    tree->SetBranchAddress("E", &E_data);

    // Create a 2D histogram
    TH2F *EX1histogram = new TH2F("EX1histogram", "E-x1", 1500, 50, 560, 400, 0.5, 16.5); // Adjust binning and range as needed

    // Fill the histogram from the tree
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        EX1histogram->Fill(x1_data, E_data);
    }

    TH2F *X1X4histogram = new TH2F("X1X4histogram", "x1-x4", 7000, 50, 480, 10000, 50, 560); // Adjust binning and range as needed

    // Fill histogram
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        X1X4histogram->Fill(x1_data, x4_data);
    }

    TH2F *X1Y4histogram = new TH2F("X1Y4histogram", "x1-y4", 2500, 50, 560, 13, 4249.5125, 4262.5125); // Adjust binning and range as needed

    // Fill histogram
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        X1Y4histogram->Fill(x1_data, y4_data);
    }

    TH2F *Z1Z4histogram = new TH2F("Z1Z4histogram", "z1-z4", 200, -6, 6, 200, -6, 6); // Adjust binning and range as needed

    // Fill histogram
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        Z1Z4histogram->Fill(z1_data, z4_data);
    }
// remove drawing for merge

    TCanvas *canvas1 = new TCanvas("canvas1", "E-x1", 800, 600);
    EX1histogram->Draw("colz");
    canvas1->Update();
    canvas1->Modified();
    canvas1->Draw();
    TCanvas *canvas2 = new TCanvas("canvas2", "x1-x4", 800, 600);
    X1X4histogram->Draw("colz");
    TGraph *q1EX1Graph = CalculatePercentileCurve(X1X4histogram, 5.);
    TGraph *q50EX1Graph = CalculatePercentileCurve(X1X4histogram, 50.);
    TGraph *q99EX1Graph = CalculatePercentileCurve(X1X4histogram, 95.);
    q1EX1Graph->SetLineColor(kRed);
    q1EX1Graph->SetLineWidth(1);
    q50EX1Graph->SetLineColor(kBlue);
    q50EX1Graph->SetLineWidth(2);
    q99EX1Graph->SetLineColor(kRed);
    q99EX1Graph->SetLineWidth(1);
    q1EX1Graph->Draw("sameL");
    q50EX1Graph->Draw("sameL");
    q99EX1Graph->Draw("sameL");
    canvas2->Update();
    canvas2->Modified();
    canvas2->Draw();
    TCanvas *canvas3 = new TCanvas("canvas3", "x1-y4", 800, 600);
    X1Y4histogram->Draw("colz");
    canvas3->Update();
    canvas3->Modified();
    canvas3->Draw();
    TCanvas *canvas4 = new TCanvas("canvas4", "z1-z4", 800, 600);
    Z1Z4histogram->Draw("colz");
    canvas4->Update();
    canvas4->Modified();
    canvas4->Draw();

    TBrowser *browser = new TBrowser("TBrowser", file);

    browser->SetTitle("ROOT Browser");
    auto EX1_lookup = makeLookUp(EX1histogram,"EX1_lookup_table.txt");
    auto X1X4_lookup = makeLookUp(X1X4histogram,"X1X4_lookup_table.txt");
    auto X1Y4_lookup = makeLookUp(X1Y4histogram,"X1Y4_lookup_table.txt",true);
    auto Z1Z4_lookup = makeLookUp(Z1Z4histogram,"Z1Z4_lookup_table.txt");

    // Clean up
//    delete browser;
//    delete canvas;
//    delete histogram;
//    delete file;
}