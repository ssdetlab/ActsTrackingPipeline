#pragma once

#include "ActsLUXEPipeline/LUXEMeasurementsCreator.hpp"
#include "ActsLUXEPipeline/LUXEGeometryConstraints.hpp"
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TBrowser.h>

void analyzeSeeds(std::vector<LUXETrackFinding::Seed> seeds, std::string filename) {
    TFile *file = new TFile(&filename[0], "RECREATE");
    TTree *l1tree = new TTree("L1Tree", "SeedParameters");
    TTree *l2tree = new TTree("L2Tree", "SeedParameters");
    TTree *l3tree = new TTree("L3Tree", "SeedParameters");
    struct ROOTSeed {
        float dx;
        float dz;
        float x1;
    };

    ROOTSeed s;

    l1tree->Branch("dx1", &s.dx);
    l2tree->Branch("dx2", &s.dx);
    l3tree->Branch("dx3", &s.dx);
    l1tree->Branch("dz1", &s.dz);
    l2tree->Branch("dz2", &s.dz);
    l3tree->Branch("dz3", &s.dz);
    l1tree->Branch("x1", &s.x1);
    l2tree->Branch("x1", &s.x1);
    l3tree->Branch("x1", &s.x1);

    for (auto &seed: seeds) {
        for (unsigned int l = 0; l < seed.originSourceLinks.size(); l++) {
            int layer = static_cast<int>(seed.originSourceLinks[l].geometryId().sensitive()/10);
            if (layer==2 || layer==3) {
                s.dx = std::abs(seed.distances[l-1][0]);
                s.dz = std::abs(seed.distances[l-1][2]);
                s.x1 = seed.x1;
                l1tree->Fill();
            } else if (layer==4 || layer==5) {
                s.dx = std::abs(seed.distances[l-1][0]);
                s.dz = std::abs(seed.distances[l-1][2]);
                s.x1 = seed.x1;
                l2tree->Fill();
            } else if (layer==6 || layer==7) {
                s.dx = std::abs(seed.distances[l-1][0]);
                s.dz = std::abs(seed.distances[l-1][2]);
                s.x1 = seed.x1;
                l3tree->Fill();
            }
        }
    }

    file->Write();
    file->Close();
    delete file;
}

void HistogramDatawriter(std::vector<LUXENavigator::Measurement> results, std::string filename, LUXEGeometry::GeometryOptions gOpt) {
    TFile *file = new TFile(&filename[0],"RECREATE");
    TTree *tree = new TTree("tree", "TruthParameters");
    struct ROOTMeasurement {
        unsigned int id;
        unsigned int layer;
        float x1;
        float z1;
        float x4;
        float y4;
        float z4;
        float E;
    };

    ROOTMeasurement s;

    tree->Branch("id", &s.id, "id/I");
    tree->Branch("layer", &s.layer, "id/I");
    tree->Branch("x1", &s.x1);
    tree->Branch("z1", &s.z1);
    tree->Branch("x4", &s.x4);
    tree->Branch("y4", &s.y4);
    tree->Branch("z4", &s.z4);
    tree->Branch("E", &s.E);

    for (auto& result:results) {
//        std::cout<<"globals size in results loop"<<result.globalPosition.size()<<std::endl;
        if (result.globalPosition.size()==0) {
            continue;
        }
        if (result.globalPosition[0][1]>=gOpt.layerZ.at(3) || (result.globalPosition[0][1]>=gOpt.layerZ.at(0) &&
                                                    result.globalPosition[0][0]>gOpt.chipTranslationXOdd.at(0))) {
            continue;
        }

        for (unsigned int l=0;l<result.truthParameters.size();l++) {
            s.id = result.eventId;
            s.layer = l+1;
            s.x1 = result.globalPosition[0][0];
            s.z1 = result.globalPosition[0][2];
            s.E = std::hypot(std::pow(result.truthParameters[l][4],-1),0.000511);
            s.x4 = result.globalPosition[result.truthParameters.size()-1][0];
            s.y4 = result.globalPosition[result.truthParameters.size()-1][1];
            s.z4 = result.globalPosition[result.truthParameters.size()-1][2];
            tree->Fill();
        }
    }

    file->Write();
    file->Close();
    delete file;
}

//for (size_t j = 0; j < sourceLinks.size(); j++) {
//
//Acts::Vector3 seedCandidatePos = SA(sourceLinks[j])->
//        localToGlobal(gctx, sourceLinks[j].parameters, yHat);
//if (seedCandidatePos[1] > gOpt.layerZPositions[1] - gOpt.deltaZ) {
//Acts::Vector3 d{ipParams[1], ipParams[2], ipParams[3]};
//Acts::Vector3 dHat = d/std::sqrt(d.dot(d));
//Acts::Vector3 diff = std::sqrt(1+std::pow((x4-x1)/(y4-y1),2)+std::pow((z4-z1)/(y4-y1),2))*
//                     (seedCandidatePos[1]-y1)*dHat;
//Acts::Vector3 AP = seedCandidatePos - globalPos - diff;
//
//// distance between point P and line defined by vector d and point A is
//// magnitude(AP x d) / magnitude(d)
//// project onto y axis
////                            Acts::Vector3 v = d.dot(yHat)*yHat;
////                            Acts::Vector3 APxv{v.cross(AP)};
////                            Acts::ActsScalar APdotV = AP.dot(v);
////                            Scalar distance =  std::sqrt(APxv.dot(APxv))/ pMagnitude;
//
//Scalar distance = std::sqrt(AP.dot(AP));
//
////                            std::cout<<"global position: "<<SA(sourceLinks[j])->
////                                    localToGlobal(gctx, sourceLinks[j].parameters, Acts::Vector3{0,1,0})<<std::endl;
//auto index = static_cast<int>((seedCandidatePos[1] - gOpt.layerZPositions[0] + gOpt.deltaZ)/100-1);
//Scalar delta = roadWidth[index];
//if (sourceLinks[j].eventId == sourceLinks[i].eventId) {
//std::cout << "comparing distance in " << i << "," << j
//<< " position " << distance << " " << delta << std::endl;
//std::cout<<"with y: "<<seedCandidatePos[1]<<std::endl;
//originSourceLinks.push_back(sourceLinks[j]);
//seedDistances.push_back(distance);
//}
//
//if (distance < delta) {
//seedSourceLinks.push_back(sourceLinks[j]);
//}
//}
//}
////                    if (seedSourceLinks.size()>1) {
//Seed addSeed{originSourceLinks, x1, seedDistances, seedSourceLinks, ipParams};
//seeds.push_back(addSeed);
////                    }
//}
