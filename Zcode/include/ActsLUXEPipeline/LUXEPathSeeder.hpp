#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Detector/Detector.hpp"
#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXEMeasurement.hpp"
#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"

#include <fstream>
#include <unordered_map>
#include <vector>
#include <TH2F.h>
#include <chrono>

/// @brief The ideal seeder for the LUXE simulation
/// takes the the SimMeasurements and converts them 
/// into seeds
namespace LUXETrackFinding {
    using namespace Acts::UnitLiterals;
    using Scalar = Acts::ActsScalar;
    struct Seed {
        std::vector<SimpleSourceLink> originSourceLinks;
        Acts::ActsScalar x1;
        std::vector<Acts::Vector3> distances;
        /// Source links related
        /// to the seed measurements
        std::vector<SimpleSourceLink> sourceLinks;
        /// IP parameters
        Acts::Vector4 ipParameters;
    };

    std::unordered_map <Scalar, Scalar> readLookup(std::string file) {
        std::ifstream lookupFile(&file[0]);
        std::unordered_map <Scalar, Scalar> lookupTable;
        Scalar x, y;
        while (lookupFile >> x >> y) {
            lookupTable[x] = y;
        }
        lookupFile.close();
        return lookupTable;
    }

    std::unordered_map <Scalar, std::vector<Scalar>> read3DLookup(std::string file) {
        std::ifstream lookupFile(&file[0]);
        std::unordered_map <Scalar, std::vector<Scalar>> lookupTable;
        Scalar x, y, z;
        while (lookupFile >> x >> y >> z) {
            std::vector<Scalar> v{y,z};
            lookupTable[x] = v;
        }
        lookupFile.close();
        return lookupTable;
    }

    using Seeds = std::vector<LUXETrackFinding::Seed>;

    template <typename T>
    T findClosestValue(std::unordered_map<Scalar, T>& lookupTable, Scalar x) {

        Scalar closestKey;
        T closestValue;
        bool first = true;

        for (const auto& entry : lookupTable) {
            if (first) {
                closestKey=entry.first;
                closestValue = entry.second;
                first = false;
                continue;
            }
            if (std::abs(entry.first - x) < std::abs(closestKey - x)) {
                closestKey = entry.first;
                closestValue = entry.second;
            }
        }
        return closestValue;
    }

    using idMap = std::unordered_map <int, std::vector<SimpleSourceLink>>;

    using str2mapMap = std::unordered_map <std::string, idMap>;
    using str2histMap = std::unordered_map <std::string, TH2F*>;

    str2mapMap binSourceLinks(const Acts::GeometryContext gctx,
                              const LUXEGeometry::GeometryOptions gOpt,
                              std::vector<SimpleSourceLink> sourceLinks,
                              SimpleSourceLink::SurfaceAccessor SA,
                              std::pair<int,int> bins,
                              str2histMap& slMap) {
        str2mapMap lookupTable;
        for (int i = 0; i<gOpt.staveZ.size()*2; i++) {
            float xMin = (i%2==0) ? gOpt.chipXEven.at(0)-gOpt.chipSizeX/2 : gOpt.chipXOdd.at(0)-gOpt.chipSizeX/2;
            float xMax = (i%2==0) ? gOpt.chipXEven.at(7)+gOpt.chipSizeX/2 : gOpt.chipXEven.at(7)+gOpt.chipSizeX/2;
            float yMin = gOpt.chipY-gOpt.chipSizeY/2;
            float yMax = gOpt.chipY+gOpt.chipSizeY/2;
            auto nAxisBinsX = bins.first;
            auto nAxisBinsY = bins.second;
            auto lName = "layer_"+std::to_string(i);
            auto mapName = lName+"_axismap";
            slMap.insert(std::make_pair(lName, new TH2F(&mapName[0], ";x[mm];y[mm];",
                                                                        nAxisBinsX,xMin,xMax,
                                                                        nAxisBinsY,yMin,yMax)));
            std::vector<SimpleSourceLink> temp;
            for(int bx=1; bx<slMap[lName]->GetNbinsX()+1; ++bx) {
                for (int by = 1; by < slMap[lName]->GetNbinsY() + 1; ++by) {
                    int bin = slMap[lName]->GetBin(bx, by);
                    auto binCenterX = slMap[lName]->GetXaxis()->GetBinCenter(bx);
                    auto binCenterY = slMap[lName]->GetYaxis()->GetBinCenter(by);
                    lookupTable[lName].insert(make_pair(bin, temp));
                }
            }
            if(i==0) i++;
        }
        for (auto& sl : sourceLinks) {
            auto id = sl.geometryId().sensitive();
            Acts::Vector3 globalPos = SA(sl)->
                    localToGlobal(gctx, sl.parameters, Acts::Vector3{0,1,0});
            int layer = static_cast<int>(id/10);
            if (layer!=1) {
                auto lName = "layer_"+std::to_string(layer);
//                std::cout<<"Filling Layer: "<<layer<<std::endl;
                double x = globalPos[0];
                double z = globalPos[2];
                int bin = slMap[lName]->FindBin(x, z);
                lookupTable[lName][bin].push_back(sl);
            }
        }
        return lookupTable;
    }

    Seeds LUXEPathSeeder(const Acts::GeometryContext gctx,
                         const LUXEGeometry::GeometryOptions gOpt,
                         std::shared_ptr<const Acts::Experimental::Detector> detector,
                         std::vector<SimpleSourceLink> sourceLinks,
                         std::string lookupDir) {

        std::cout<<"Preparing Look Up tables"<<std::endl;
        auto start = std::chrono::steady_clock::now();
        std::unordered_map <Scalar, Scalar> EX1LookUp = readLookup(lookupDir + "/EX1_lookup_table.txt");
        std::unordered_map <Scalar, Scalar> X1X4LookUp = readLookup(lookupDir + "/X1X4_lookup_table.txt");
        std::unordered_map <Scalar, Scalar> X1Y4LookUp = readLookup(lookupDir + "/X1Y4_lookup_table.txt");
        std::unordered_map <Scalar, Scalar> Z1Z4LookUp = readLookup(lookupDir + "/Z1Z4_lookup_table.txt");
        auto lookupEnd = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(lookupEnd - start);

        std::cout << "Time taken: " << duration.count() << "ms" << std::endl;

        SimpleSourceLink::SurfaceAccessor SA{*detector};

        std::pair<int, int> bins = std::make_pair(3500,200);
        str2histMap slMap;
        std::cout << "Preparing bins"<< std::endl;
        auto lookupTable = binSourceLinks(gctx, gOpt, sourceLinks, SA, bins, slMap);
        auto binningEnd = std::chrono::steady_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(binningEnd - lookupEnd);
        std::cout << "Time taken: " << duration.count() << "ms" << std::endl;
        Seeds seeds;
        Acts::Vector4 ipParams;

//        Acts::Vector4 roadWidthX{100_um,250_um,350_um,450_um};
//        Acts::Vector4 roadWidthZ{50_um,80_um,120_um,150_um};

        Acts::Vector4 roadWidthX{1200_um,2000_um,2000_um,2000_um};
        Acts::Vector4 roadWidthZ{1200_um,2000_um,2000_um,2000_um};

        start = std::chrono::steady_clock::now();
        auto end = std::chrono::steady_clock::now();
        for (size_t i=0; i<sourceLinks.size(); i++) {
            Acts::Vector3 globalPos = SA(sourceLinks[i])->
                localToGlobal(gctx, sourceLinks[i].parameters, Acts::Vector3{0,1,0});

            Scalar x1 = globalPos[0];
            Scalar y1 = globalPos[1];
            Scalar z1 = globalPos[2];

            if (i%(sourceLinks.size()/10)==0) {
                end = std::chrono::steady_clock::now();
                duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
                std::cout << "Completed "<<(i*100)/sourceLinks.size()<<"%"<< std::endl;
                std::cout << "Total time elapsed: " << duration.count()/1000 << "s" << std::endl;
            }
            if (y1 <= gOpt.layerZPositions.at(0) ||
                    (y1 < gOpt.layerZPositions.at(1)-gOpt.deltaZ &&
                     x1 <= gOpt.chipXOdd.at(0)-gOpt.chipSizeX/2+1.8)) { // I can explain the 1.8

                std::vector<SimpleSourceLink> seedSourceLinks{sourceLinks[i]};
                std::vector<Acts::Vector3> seedDistances;
                std::vector<SimpleSourceLink> originSourceLinks{sourceLinks[i]};

                for (size_t j = i + 1; j < sourceLinks.size(); j++) {
                    if (sourceLinks[i].eventId == sourceLinks[j].eventId) {
                        originSourceLinks.push_back(sourceLinks[j]);
                    }
                }
                if (originSourceLinks.size() > 3) {

                    Scalar E = LUXETrackFinding::findClosestValue(EX1LookUp, x1);
                    Scalar x4 = LUXETrackFinding::findClosestValue(X1X4LookUp, x1);
                    Scalar y4 = LUXETrackFinding::findClosestValue(X1Y4LookUp, x1);
                    y4 = (y4 < gOpt.layerZPositions.at(3)) ? gOpt.staveZ.at(7) : gOpt.staveZ.at(6); // get rid of lookup binning errors
                    Scalar z4 = LUXETrackFinding::findClosestValue(Z1Z4LookUp, z1);
                    Scalar electron_mass = 0.000511;
                    Scalar pMagnitude = std::sqrt(std::pow(E, 2) -
                                                  std::pow(electron_mass, 2));
                    Scalar vMagnitude = pMagnitude / std::sqrt(std::pow((x4 - x1), 2) +
                                                               std::pow((y4 - y1), 2) +
                                                               std::pow((z4 - z1), 2));
                    ipParams = {E,
                                vMagnitude * (x4 - x1),
                                vMagnitude * (y4 - y1),
                                vMagnitude * (z4 - z1)};

                    Acts::Vector3 d{ipParams[1], ipParams[2], ipParams[3]};
                    Acts::Vector3 dHat = d / std::sqrt(d.dot(d));

//                    roadWidthX[1] = (0.3*x1+85) * 1_um;
//                    roadWidthX[2] = (0.8*x1+140) * 1_um;
//                    roadWidthX[3] = (0.9*x1+195) * 1_um;

                    std::vector<Acts::Vector3> layerPointers;
                    for (int il = 0; il < gOpt.staveZ.size(); il++) {
                            Acts::Vector3 Lpos = std::sqrt(1 + std::pow((x4 - x1) / (y4 - y1), 2) +
                                                           std::pow((z4 - z1) / (y4 - y1), 2)) *
                                                 (gOpt.staveZ.at(il)-y1) * dHat;
                            layerPointers.push_back(Lpos);
                    }

                    for (int k = 0; k < layerPointers.size(); k++) {
                        if (k==1) continue;
                        if (k==0 && (x1 > gOpt.chipXEven.at(8) + gOpt.chipSizeX/2 ||
                                     x1<=gOpt.chipXOdd.at(0))-gOpt.chipSizeX/2 + 1.8) continue;

                        Acts::Vector3 refPoint = globalPos + layerPointers[k];

                        std::pair<float, float> Xlim = (k % 2 == 0) ?
                                                       std::make_pair(gOpt.chipXEven.at(0) - gOpt.chipSizeX/2,
                                                                      gOpt.chipXEven.at(8) + gOpt.chipSizeX/2) :
                                                       std::make_pair(gOpt.chipXOdd.at(0) - gOpt.chipSizeX/2,
                                                                      gOpt.chipXOdd.at(8) + gOpt.chipSizeX/2);
                        std::pair<float, float> Zlim = std::make_pair(gOpt.chipY,
                                                                      gOpt.chipY + gOpt.chipSizeY);

                        Scalar topX = std::min(refPoint[0] + roadWidthX[k / 2],
                                               static_cast<double>(Xlim.second));
                        Scalar botX = std::max(refPoint[0] - roadWidthX[k / 2],
                                               static_cast<double>(Xlim.first));
                        Scalar topZ = std::min(refPoint[2] + roadWidthZ[k / 2],
                                               static_cast<double>(Zlim.second));
                        Scalar botZ = std::max(refPoint[2] - roadWidthZ[k / 2],
                                               static_cast<double>(Zlim.first));

                        std::string lName = "layer_" + std::to_string(k);

                        //========FOR TESTING========
                        for (auto og : originSourceLinks) {
                            unsigned int gId = og.geometryId().sensitive()/10;
                            if (gId == 1) {
                                continue;
                            }
                            if (k==gId) {
                                Acts::Vector3 ogPos = SA(og) -> localToGlobal(gctx, og.parameters, Acts::Vector3{0,1,0});
                                seedDistances.push_back(ogPos-refPoint);
                            }
                        }
                        //===========================
                        int botLeftBin =  slMap[&lName[0]]->FindBin(botX, botZ-0.001);
                        int rightBin =    slMap[&lName[0]]->FindBin(topX, botZ-0.001);
                        int topRightBin = slMap[&lName[0]]->FindBin(topX, topZ+0.001);

                        int pathWidthInBins = rightBin - botLeftBin;
                        int currentBin = botLeftBin;
                        while (currentBin <= topRightBin + 1) {
                            while (currentBin <= rightBin + 1) {
                                auto sourceLinksToAdd = lookupTable[&lName[0]][currentBin];
                                seedSourceLinks.insert(seedSourceLinks.end(),
                                                       sourceLinksToAdd.begin(),
                                                       sourceLinksToAdd.end());
                                currentBin++;
                            }
                            currentBin += bins.first + 2 - (pathWidthInBins + 1); // explicit
                            rightBin += bins.first + 2;
                        }
                    }
                    Seed seed{originSourceLinks, x1, seedDistances, seedSourceLinks, ipParams};
                    seeds.push_back(seed);
                }
            }
        }
        return seeds;
    }
} // LUXETrackFinding

