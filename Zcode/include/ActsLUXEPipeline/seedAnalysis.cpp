#pragma once

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TBrowser.h>
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

//std::unordered_map<float, float> makeLookUp(TH2F* histogram, std::string filename) {
//
//    std::unordered_map<float, float> lookupTable;
//    // Loop over bins to build the lookup table
//    for (int binx = 1; binx <= histogram->GetNbinsX(); ++binx) {
//        float x = histogram->GetXaxis()->GetBinCenter(binx);
//        float cumSum = 0.0;
//        std::vector<float> values;
//        for (int biny = 1; biny <= histogram->GetNbinsY(); ++biny) {
//            cumSum += histogram->GetBinContent(binx, biny);
//            values.push_back(cumSum);
//        }
//        float totalSum = cumSum;
//        auto it = std::lower_bound(values.begin(), values.end(), cumSum*0.7);
//
//        float E_median = histogram->GetYaxis()->GetBinCenter(std::distance(values.begin(), it));
//        lookupTable[x] = E_median;
//    }
//    FILE *lookupFile = fopen(&filename[0], "w");
//    if (lookupFile) {
//        for (const auto &entry : lookupTable) {
//            fprintf(lookupFile, "%f %f\n", entry.first, entry.second);
//        }
//        fclose(lookupFile);
//    } else {
//        std::cerr << "Failed to open lookup table file for writing." << std::endl;
//    }
//    return lookupTable;
//}

void seedAnalysis() {
    TFile *file = new TFile("seed_data.root");

    // Retrieve the tree from the file
    TTree *l1tree;
    TTree *l2tree;
    TTree *l3tree;
    file->GetObject("L1Tree", l1tree);
    file->GetObject("L2Tree", l2tree);
    file->GetObject("L3Tree", l3tree);

    Float_t dx1_data, dz1_data, dx2_data, dz2_data, dx3_data, dz3_data, x11_data, x12_data, x13_data;
    l1tree->SetBranchAddress("dx1", &dx1_data);
    l1tree->SetBranchAddress("dz1", &dz1_data);
    l1tree->SetBranchAddress("x1", &x11_data);

    l2tree->SetBranchAddress("dx2", &dx2_data);
    l2tree->SetBranchAddress("dz2", &dz2_data);
    l2tree->SetBranchAddress("x1", &x12_data);

    l3tree->SetBranchAddress("dx3", &dx3_data);
    l3tree->SetBranchAddress("dz3", &dz3_data);
    l3tree->SetBranchAddress("x1", &x13_data);

    // Create a 2D histogram
    TH2F *X1DX1histogram = new TH2F("X1DX1histogram", "x1-dx1", 400, 50, 550,200,0,1.5); // Adjust binning and range as needed

    // Fill the histogram from the tree
    Long64_t nEntries = l1tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        l1tree->GetEntry(i);
        X1DX1histogram->Fill(x11_data, dx1_data);
    }

    TH2F *X1DZ1histogram = new TH2F("X1DZ1histogram", "x1-dz1", 400, 50, 550,200,0,1.5); // Adjust binning and range as needed

    // Fill the histogram from the tree
    for (Long64_t i = 0; i < nEntries; ++i) {
        l1tree->GetEntry(i);
        X1DZ1histogram->Fill(x11_data, dz1_data);
    }

    TH2F *X1DX2histogram = new TH2F("X1DX2histogram", "x1-dx2", 400, 50, 550,200,0,1.5); // Adjust binning and range as needed

    // Fill histogram
    nEntries = l2tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        l2tree->GetEntry(i);
        X1DX2histogram->Fill(x12_data, dx2_data);
    }

    TH2F *X1DZ2histogram = new TH2F("X1DZ2histogram", "x1-dz2", 400, 50, 550,200,0,1.5); // Adjust binning and range as needed

    // Fill histogram
    for (Long64_t i = 0; i < nEntries; ++i) {
        l2tree->GetEntry(i);
        X1DZ2histogram->Fill(x12_data, dz2_data);
    }

    TH2F *X1DX3histogram = new TH2F("X1DX3histogram", "x1-dx3", 400, 50, 550,200,0,1.5); // Adjust binning and range as needed

    // Fill histogram
    nEntries = l3tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        l3tree->GetEntry(i);
        X1DX3histogram->Fill(x13_data, dx3_data);
    }

    TH2F *X1DZ3histogram = new TH2F("X1DZ3histogram", "x1-dz3", 400, 50, 550,200,0,2); // Adjust binning and range as needed

    // Fill histogram
    for (Long64_t i = 0; i < nEntries; ++i) {
        l3tree->GetEntry(i);
        X1DZ3histogram->Fill(x13_data, dz3_data);
    }

    TCanvas *l1_x = new TCanvas("l1_x", "x1-dx1", 800, 600);
    X1DX1histogram->Draw("colz");
    l1_x->Update();
    l1_x->Modified();
    l1_x->Draw();
    TGraph *q95X1DX1Graph = CalculatePercentileCurve(X1DX1histogram, 95.);
    q95X1DX1Graph->SetLineColor(kRed);
    q95X1DX1Graph->SetLineWidth(2);
    q95X1DX1Graph->Draw("sameL");

    TCanvas *l1_z = new TCanvas("l1_z", "x1-dz1", 800, 600);
    X1DZ1histogram->Draw("colz");
    l1_z->Update();
    l1_z->Modified();
    l1_z->Draw();
    TGraph *q95X1DZ1Graph = CalculatePercentileCurve(X1DZ1histogram, 95.);
    q95X1DZ1Graph->SetLineColor(kRed);
    q95X1DZ1Graph->SetLineWidth(2);
    q95X1DZ1Graph->Draw("sameL");

    TCanvas *l2_x = new TCanvas("l2_x", "x1-dx2", 800, 600);
    X1DX2histogram->Draw("colz");
    l2_x->Update();
    l2_x->Modified();
    l2_x->Draw();
    TGraph *q95X1DX2Graph = CalculatePercentileCurve(X1DX2histogram, 95.);
    q95X1DX2Graph->SetLineColor(kRed);
    q95X1DX2Graph->SetLineWidth(2);
    q95X1DX2Graph->Draw("sameL");

    TCanvas *l2_z = new TCanvas("l2_z", "x1-dz2", 800, 600);
    X1DZ2histogram->Draw("colz");
    l2_z->Update();
    l2_z->Modified();
    l2_z->Draw();
    TGraph *q95X1DZ2Graph = CalculatePercentileCurve(X1DZ2histogram, 95.);
    q95X1DZ2Graph->SetLineColor(kRed);
    q95X1DZ2Graph->SetLineWidth(2);
    q95X1DZ2Graph->Draw("sameL");

    TCanvas *l3_x = new TCanvas("l3_x", "x1-dx3", 800, 600);
    X1DX3histogram->Draw("colz");
    l3_x->Update();
    l3_x->Modified();
    l3_x->Draw();
    TGraph *q95X1DX3Graph = CalculatePercentileCurve(X1DX3histogram, 95);
    q95X1DX3Graph->SetLineColor(kRed);
    q95X1DX3Graph->SetLineWidth(2);
    q95X1DX3Graph->Draw("sameL");

    TCanvas *l3_z = new TCanvas("l3_z", "x1-dz3", 800, 600);
    X1DZ3histogram->Draw("colz");
    l3_z->Update();
    l3_z->Modified();
    l3_z->Draw();
    TGraph *q95X1DZ3Graph = CalculatePercentileCurve(X1DZ3histogram, 95);
    q95X1DZ3Graph->SetLineColor(kRed);
    q95X1DZ3Graph->SetLineWidth(2);
    q95X1DZ3Graph->Draw("sameL");

    TBrowser *browser = new TBrowser("TBrowser", file);

    browser->SetTitle("ROOT Browser");
//    auto EX1_lookup = makeLookUp(EX1histogram,"ED1_lookup_table.txt");
//    auto X1X4_lookup = makeLookUp(X1X4histogram,"ED2_lookup_table.txt");
//    auto Z1Z4_lookup = makeLookUp(Z1Z4histogram,"ED3_lookup_table.txt");

    // Clean up
//    delete browser;
//    delete canvas;
//    delete histogram;
//    delete file;
}