#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Detector/Detector.hpp"
#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"
//#include "ActsLUXEPipeline/LUXEDataContainers.hpp"

#include <fstream>
#include <unordered_map>

/// @brief The ideal seeder for the LUXE simulation
/// takes the the SimMeasurements and converts them 
/// into seeds
namespace LUXETrackFinding {
    using namespace Acts::UnitLiterals;
    using Scalar = Acts::ActsScalar;
    struct Seed {
        SimpleSourceLink originSourceLink;
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
        return closestValue;
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
        Seeds seeds;
        Acts::Vector4 ipParams;
        Scalar delta = 500_um;
        for (size_t i=0; i<sourceLinks.size(); i++) {
            Acts::Vector3 globalPos = SA(sourceLinks[i])->
                localToGlobal(gctx, sourceLinks[i].parameters, Acts::Vector3{0,1,0});
            // TODO: look for an elegant way to realise if the sourceLink is in L1
            if (globalPos[1] <= gOpt.layerZPositions[1]-gOpt.deltaZ) {
                Scalar x1 = globalPos[0];
                Scalar y1 = globalPos[1];
                Scalar z1 = 1e13*globalPos[2]; //
                Scalar E = LUXETrackFinding::findClosestValue(EX1LookUp,x1);
                Scalar x4 = LUXETrackFinding::findClosestValue(X1X4LookUp,x1);
                Scalar y4 = gOpt.layerZPositions[3]-gOpt.deltaZ;
                Scalar z4 = LUXETrackFinding::findClosestValue(Z1Z4LookUp,z1);
                Scalar electron_mass = 0.000511;
                z1 = 1e-13*z1;
                Scalar pMagnitude = std::sqrt(std::pow(E,2)-
                                              std::pow(electron_mass,2));
                Scalar vMagnitude = pMagnitude/std::sqrt(std::pow((x4-x1),2)+
                                             std::pow((y4-y1),2)+
                                             std::pow((z4-z1),2));
                ipParams = {E,
                          vMagnitude*(x4-x1),
                          vMagnitude*(y4-y1),
                          vMagnitude*(z4-z1)};

                std::cout<<ipParams<<std::endl;
                std::vector<SimpleSourceLink> seedSourceLinks;

                for (size_t j=0; j<sourceLinks.size(); j++) {
                    //TODO: EFFICIENCY
                    Acts::Vector3 seedCandidatePos = SA(sourceLinks[j])->
                            localToGlobal(gctx, sourceLinks[j].parameters, Acts::Vector3{0,1,0});
                    if (seedCandidatePos[1] > gOpt.layerZPositions[1]-gOpt.deltaZ) {
                        Acts::Vector3 diff = seedCandidatePos-globalPos;
                        // distance between point P and line defined by vector d and point A is
                        // magnitude(AP x d) / magnitude(d)
                        Acts::Vector3 d{ipParams[1],ipParams[2],ipParams[3]};
                        Acts::Vector3 APxd{d.cross(diff)};
//                            ipParams[2]*diff[0]-ipParams[1]*diff[1],
//                                                   ipParams[1]*diff[2]-ipParams[3]*diff[0],
//                                                   ipParams[3]*diff[1]-ipParams[2]*diff[2]};
                        Scalar distance = std::sqrt(std::pow(APxd[0],2)+
                                                    std::pow(APxd[1],2)+
                                                    std::pow(APxd[2],2))/pMagnitude;
                        std::cout<<"comparing distance in "<<i<<","<<j
                        <<" position "<<distance<<" "<<delta<<std::endl;
                        if (distance<delta) {
                            seedSourceLinks.push_back(sourceLinks[j]);
                        }
                    }
                }
                Seed addSeed{sourceLinks[i],seedSourceLinks,ipParams};
                seeds.push_back(addSeed);
            }
        }
        return seeds;
    }
}

