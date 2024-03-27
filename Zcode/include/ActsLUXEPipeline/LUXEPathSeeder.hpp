#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Detector/Detector.hpp"
#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/LUXEGeometry.hpp"

#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"
//#include "ActsLUXEPipeline/LUXEDataContainers.hpp"

//#include "ActsLUXEPipeline/LUXENavigator.hpp"

#include <fstream>
#include <unordered_map>
#include <TH2F.h>

/// @brief The ideal seeder for the LUXE simulation
/// takes the the SimMeasurements and converts them 
/// into seeds
namespace LUXETrackFinding {
    using namespace Acts::UnitLiterals;
    using Scalar = Acts::ActsScalar;
    struct Seed {
        std::vector<SimpleSourceLink> originSourceLinks;
        Acts::ActsScalar x1;
        std::vector<Acts::ActsScalar> distances;
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

    using Seeds = std::vector<LUXETrackFinding::Seed>;

    Scalar findClosestValue(std::unordered_map<Scalar, Scalar>& lookupTable, Scalar x) {

        Scalar closestKey;
        Scalar closestValue;
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
        if (closestValue>40) {
            std::cout<<"Closest value was: "<<closestKey<<std::endl;
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
        for (int i = 2; i<gOpt.layerZPositions.size()*2; i++) {
            float xMin = (i%2==0) ? gOpt.chipTranslationXEven.at(0) : gOpt.chipTranslationXOdd.at(0);
            float xMax = (i%2==0) ? gOpt.chipTranslationXEven.at(8)+gOpt.chipSizeX : gOpt.chipTranslationXOdd.at(8)+gOpt.chipSizeX;
            float yMin = gOpt.chipTranslationY;
            float yMax = gOpt.chipTranslationY+gOpt.chipSizeY;
            int nAxisBinsX = bins.first;
            int nAxisBinsY = bins.second;
            auto lName = "layer_"+std::to_string(i);
            std::cout<<lName<<"LAYERCHECK"<<std::endl;
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
        }
        for (auto& sl : sourceLinks) {
            auto id = sl.geometryId().sensitive();
            Acts::Vector3 globalPos = SA(sl)->
                    localToGlobal(gctx, sl.parameters, Acts::Vector3{0,1,0});
            int layer = static_cast<int>(id/10);
            if (layer>=2) {
                auto lName = "layer_"+std::to_string(layer);
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
        std::unordered_map <Scalar, Scalar> EX1LookUp = readLookup(lookupDir + "/EX1_lookup_table.txt");
        std::unordered_map <Scalar, Scalar> X1X4LookUp = readLookup(lookupDir + "/X1X4_lookup_table.txt");
        std::unordered_map <Scalar, Scalar> Z1Z4LookUp = readLookup(lookupDir + "/Z1Z4_lookup_table.txt");

        SimpleSourceLink::SurfaceAccessor SA{*detector};

        std::pair<int, int> bins = std::make_pair(30,10);
        str2histMap slMap;
        auto lookupTable = binSourceLinks(gctx, gOpt, sourceLinks, SA, bins, slMap);

        Seeds seeds;
        Acts::Vector4 ipParams;
        Acts::Vector3 roadWidthX{1000_um,1500_um,2500_um};
        Acts::Vector3 roadWidthZ{1500_um,1500_um,1550_um};

        for (size_t i=0; i<sourceLinks.size(); i++) {
            Acts::Vector3 globalPos = SA(sourceLinks[i])->
                localToGlobal(gctx, sourceLinks[i].parameters, Acts::Vector3{0,1,0});
            // TODO: look for an elegant way to realise if the sourceLink is in L1
            if (globalPos[1] <= gOpt.layerZPositions[0] ||
                    (globalPos[1] <= gOpt.layerZPositions[1]-gOpt.deltaZ &&
                     globalPos[0] <= gOpt.chipTranslationXOdd.at(0))) {
                Scalar x1 = globalPos[0];
                Scalar y1 = globalPos[1];
                Scalar z1 = globalPos[2];
                Scalar E = LUXETrackFinding::findClosestValue(EX1LookUp,x1);
                Scalar x4 = LUXETrackFinding::findClosestValue(X1X4LookUp,x1);
                Scalar y4 = (x4>gOpt.chipTranslationXEven.at(8)+gOpt.chipSizeX) ?
                                          gOpt.layerZPositions[3]-gOpt.deltaZ :
                                          gOpt.layerZPositions[3]+gOpt.deltaZ;

                Scalar z4 = LUXETrackFinding::findClosestValue(Z1Z4LookUp,z1);
                Scalar electron_mass = 0.000511;
                Scalar pMagnitude = std::sqrt(std::pow(E,2)-
                                              std::pow(electron_mass,2));
                Acts::Vector3 yHat{0,1,0};

                std::vector<SimpleSourceLink> seedSourceLinks{sourceLinks[i]};
                std::vector<Acts::ActsScalar> seedDistances;
                std::vector<SimpleSourceLink> originSourceLinks{sourceLinks[i]};

                Scalar vMagnitude = pMagnitude/std::sqrt(std::pow((x4-x1),2)+
                                                         std::pow((y4-y1),2)+
                                                         std::pow((z4-z1),2));
                ipParams = {E,
                            vMagnitude*(x4-x1),
                            vMagnitude*(y4-y1),
                            vMagnitude*(z4-z1)};
                Acts::Vector3 d{ipParams[1], ipParams[2], ipParams[3]};
                Acts::Vector3 dHat = d/std::sqrt(d.dot(d));

                std::vector<double> factor{0,2};
                Scalar layerGap = gOpt.layerZPositions[1]-gOpt.layerZPositions[0];
                std::vector<Acts::Vector3> layerPointers;
                for (int il=1; il<gOpt.layerZPositions.size();il++) {
                    for (auto f : factor) {
                        Acts::Vector3 Lpos = (y1<gOpt.layerZPositions[0]) ?
                                                std::sqrt(1+std::pow((x4-x1)/(y4-y1),2)+std::pow((z4-z1)/(y4-y1),2))*
                                                (il*layerGap + (2-f)*gOpt.deltaZ)*dHat : // TODO: maybe static cast i
                                                std::sqrt(1+std::pow((x4-x1)/(y4-y1),2)+std::pow((z4-z1)/(y4-y1),2))*
                                                (il*layerGap - f*gOpt.deltaZ)*dHat ;
                        layerPointers.push_back(Lpos);
                    }
                }
                for (int k=0 ; k<layerPointers.size(); k++) {

                    Acts::Vector3 refPoint = globalPos + layerPointers[k];
                    std::pair<float,float> Xlim = (k%2==0) ?
                            std::make_pair(gOpt.chipTranslationXEven.at(0),
                                           gOpt.chipTranslationXEven.at(8)+gOpt.chipSizeX) :
                            std::make_pair(gOpt.chipTranslationXOdd.at(0),
                                           gOpt.chipTranslationXOdd.at(8)+gOpt.chipSizeX);
                    std::pair<float,float> Zlim = std::make_pair(gOpt.chipTranslationY,
                                                                 gOpt.chipTranslationY+gOpt.chipSizeY);
                    Scalar topX = std::min(refPoint[0] + roadWidthX[k/2],
                                           static_cast<double>(Xlim.second));
                    Scalar botX = std::max(refPoint[0] - roadWidthX[k/2],
                                           static_cast<double>(Xlim.first));
                    Scalar topZ = std::min(refPoint[2] + roadWidthZ[k/2],
                                           static_cast<double>(Zlim.second));
                    Scalar botZ = std::max(refPoint[2] - roadWidthZ[k/2],
                                           static_cast<double>(Zlim.first));
                    if (botX>=topX || botZ>=topZ) {
//                            std::cout<<"refPoint: "<<refPoint<<std::endl;
//                            std::cout<<"Xlim.first: "<<Xlim.first<<std::endl;
//                            std::cout<<"Xlim.second: "<<Xlim.second<<std::endl;
//                            std::cout<<"Zlim.first: "<<Zlim.first<<std::endl;
//                            std::cout<<"Zlim.second: "<<Zlim.second<<std::endl;
//                            std::cout<<"4 corners: "<<topX<<" "<<botX<<" "<<topZ<<" "<<botZ<<std::endl;
                        std::cout<<"path doesn't contain detector layer at pointer "<<i<<std::endl;
                    }
                    std::string lName = "layer_"+std::to_string(k+2);

                    int botLeftBin  = slMap[&lName[0]]->FindBin(botX,botZ);
                    int rightBin    = slMap[&lName[0]]->FindBin(topX,botZ);
                    int topRightBin = slMap[&lName[0]]->FindBin(topX,topZ);

                    int pathWidthInBins = rightBin-botLeftBin;
                    int currentBin = botLeftBin;
                    while (currentBin <= topRightBin) {
                        while (currentBin <= rightBin) {
                            auto sourceLinksToAdd = lookupTable[&lName[0]][currentBin];
                            seedSourceLinks.insert(seedSourceLinks.end(),
                                                   sourceLinksToAdd.begin(),
                                                   sourceLinksToAdd.end());
                            currentBin++;
                        }
                        currentBin += bins.first + 2 - (pathWidthInBins+1); // explicit
                        rightBin += bins.first + 2;
                    }
                }
                for (size_t j=i+1; j<sourceLinks.size(); j++) {
                    if (sourceLinks[i].eventId == sourceLinks[j].eventId) {
                        originSourceLinks.push_back(sourceLinks[j]);
                    }
                }
                Seed seed{originSourceLinks,0.0,seedDistances,seedSourceLinks,ipParams};
                seeds.push_back(seed);
            }
        }
        return seeds;
    }
} // LUXETrackFinding

