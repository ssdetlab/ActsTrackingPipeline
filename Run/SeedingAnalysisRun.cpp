#include "ActsLUXEPipeline/Sequencer.hpp"

#include "ActsLUXEPipeline/LUXEBinnedMagneticField.hpp"
#include "ActsLUXEPipeline/ConstantBoundedField.hpp"
#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXEPathSeeder.hpp"
#include "ActsLUXEPipeline/LUXEGeometryConstraints.hpp"
#include "ActsLUXEPipeline/LUXEMagneticField.hpp"
#include "ActsLUXEPipeline/LUXENavigator.hpp"
#include "ActsLUXEPipeline/LUXEMeasurementsCreator.hpp"
#include "ActsLUXEPipeline/LUXEPathSeeder.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"

#include <filesystem>
#include <string>
#include <iostream>
#include <random>
#include <cstdlib>


#include "ActsLUXEPipeline/Utils.hpp"

/// @brief Run the propagation through
/// a uniform energy spectrum and record the
/// energy vs position histograms for each layer
int main(int argc, char* argv[]) {
    using namespace LUXENavigator;
    Acts::Logging::Level logLevel = Acts::Logging::INFO;

    // setup the sequencer first w/ config derived from options
    Sequencer::Config seqCfg;
    seqCfg.events = 10;
    seqCfg.numThreads = -1;
    Sequencer sequencer(seqCfg);

    const std::vector<std::pair<Acts::ActsScalar,Acts::ActsScalar>> MagneticFieldBounds =
            {std::make_pair(-1000_mm,1000_mm),
             std::make_pair(1450_mm,2650_mm),
             std::make_pair(-100_mm,100_mm)};

    auto transformPos = [&](const Acts::Vector3& pos) {
        for (int i=0;i<3;i++) {
            if (pos[i]<= MagneticFieldBounds[i].first ||
                pos[i]>  MagneticFieldBounds[i].second) {
                return Acts::Vector3{0,1400,0};
            }
        }
        return pos;
    };

    // map (Bx,By,Bz) -> (Bx,By,Bz)
    auto transformBField = [](const Acts::Vector3& field, const Acts::Vector3&) {
        return field;
    };
    LUXEMagneticField::vGridOptions gridOpt;
    gridOpt.xBins = {-1000,-1, 0.,200, 1000.};
    gridOpt.yBins = {1300,1400,1450,1451, 2050.,2649,2650.,2651};
    gridOpt.zBins = {-100,-99, 0.,1, 100.};

    Acts::GeometryContext gctx;
    Acts::MagneticFieldContext mctx;
    LUXEGeometry::GeometryOptions gOpt;
    double B_z = 0.95_T;

    Acts::Extent dipoleExtent;
    dipoleExtent.set(Acts::binX, -1000_mm, 1000_mm);
    dipoleExtent.set(Acts::binY, 1450_mm, 2650_mm);
    dipoleExtent.set(Acts::binZ, -100_mm, 100_mm);

    auto BField = LUXEMagneticField::buildBinnedBField(
            LUXEMagneticField::ConstantBoundedField(Acts::Vector3(0., 0., B_z), dipoleExtent),
            transformPos, transformBField, gridOpt, mctx);
    auto BFieldPtr = std::make_shared<Acts::InterpolatedBFieldMap<LUXEMagneticField::vGrid>>(BField);

    // Build the LUXE detector
    std::string gdmlPath = "lxgeomdump_ip_tracker_positron.gdml";

    std::vector<std::string> staves = {"OPPPSensitive"};
    std::vector<std::string> chamber = {"VCWindowPanel"};
    auto positronArmBpr = LUXEGeometry::makeBlueprintLUXE(gdmlPath, staves, gOpt);
    auto detector = LUXEGeometry::buildLUXEDetector(std::move(positronArmBpr), gctx, gOpt);

    Acts::ViewConfig pConfig = Acts::s_viewSensitive;
    Acts::ObjVisualization3D volumeObj;
    for (auto& vol : detector->rootVolumes()) {
        for (auto& surf : vol->surfaces()) {
                Acts::GeometryView3D::drawSurface(
                        volumeObj, *(surf), gctx,
                        Acts::Transform3::Identity(), pConfig);
            std::cout<<"Surface x transform: "<<surf->center(gctx)[0]<<std::endl;
            std::cout<<"Surface y transform: "<<surf->center(gctx)[1]<<std::endl;
            std::cout<<"Surface z transform: "<<surf->center(gctx)[2]<<std::endl;
            std::cout<<"Surface bounds: "<<surf->bounds()<<std::endl;
        }
    }

    Acts::ActsScalar m_e = 0.000511;
    std::vector<LUXENavigator::Measurement> results = LUXENavigator::loadMeasurementsFromFile("10k_measurements.dat"); //Zmeasurements/ AND EXTRA ZERO
    SimpleSourceLink::SurfaceAccessor SA{*detector};
    std::vector<SimpleSourceLink> sl4Seeding;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> fake_x(-gOpt.chipSizeX/2,gOpt.chipSizeX/2);
    std::uniform_real_distribution<> fake_y(-gOpt.chipSizeY/2,gOpt.chipSizeY/2);
    
    for (auto& result:results) {
        for (unsigned int l=0;l<result.sourceLinks.size();l++) {
            sl4Seeding.push_back(result.sourceLinks[l]);
        }
    }
    int fakesPerLayer = 0;
    Acts::SquareMatrix2 fakeCov = Acts::SquareMatrix2::Identity();

    for (auto& vol : detector->rootVolumes()) {
        for (auto& surf : vol->surfaces()) {
            for (int i=1 ; i<=fakesPerLayer ; i++) {
                SimpleSourceLink fakeSl({fake_x(gen),fake_y(gen)} , fakeCov, surf -> geometryId(), -i); //sl4Seeding.size()+
                sl4Seeding.push_back(fakeSl);
            }
        }
    }

    int fakes = fakesPerLayer*9*8;
    std::cout<<fakes<<" fakes"<<std::endl;

    int multiplicity;
    if (argc>1) {
        multiplicity = std::atoi(argv[1]);
    } else {
        multiplicity = 0;
    }
    //10,100,,10000
    auto start = sl4Seeding.begin();
    auto end = sl4Seeding.end();
    if (multiplicity) {
        end = start + multiplicity;
    }
    std::vector<LUXETrackFinding::Seed> badSeeds;
    std::vector<SimpleSourceLink> slBatch(start,end);
    std::vector<LUXETrackFinding::Seed> seeds = LUXETrackFinding::LUXEPathSeeder(gctx, gOpt, detector, slBatch,
                                                     "/Users/alonlevi/CLionProjects/LUXEPipeline/build/Zlookups/"); //

            SimpleSourceLink::SurfaceAccessor SAseed{*detector};
            int count = 1;
            pConfig.color = {0, 250, 0};
            int layerZero = 0;
            double efficiency = 0.;
            double fakeEfficiency = 0.;
            double logFakeEfficiency = 0.;
            size_t counter;
            std::vector<int> fakeCounter;
            double totalCombinations;
            bool L0;
            int index;
            int OGindex;
            for (auto seed: seeds) {
                counter = 0;
                totalCombinations = 1.;
                fakeCounter = {0, 0, 0, 0, 0, 0, 0};
                L0 = (seed.originSourceLinks[0].geometryId().sensitive() / 10 == 1);
                if (L0 || seed.originSourceLinks[0].geometryId().sensitive() / 10 == 0) {
                    layerZero++;
                }
                for (auto sl: seed.sourceLinks) {
                    if (sl.eventId > 0) {
                        index = sl.geometryId().sensitive() / 10 - 1;
                        if (L0 || index != 0) {
                            fakeCounter[index]++;
                        }
//                        if (std::count(seed.originSourceLinks.begin(), seed.originSourceLinks.end(), sl)) {
//                            counter++;
//                        }
                    }
                }
                for (auto sl: seed.originSourceLinks) {
                    if (std::count(seed.sourceLinks.begin(), seed.sourceLinks.end(), sl)) {
                        counter++;
                    }
                }
                if (counter == seed.originSourceLinks.size()) {
                    efficiency++;
                    for (auto a: fakeCounter) {
                        if (a != 0) totalCombinations = totalCombinations * static_cast<double>(a);
                    }
                    fakeEfficiency += 1 / totalCombinations;
                    logFakeEfficiency += 1 / (std::log(std::exp(1) - 1 + totalCombinations));
                } else {
                    badSeeds.push_back(seed);
                }
            }

            fakeEfficiency = fakeEfficiency / efficiency * 100;
            logFakeEfficiency = logFakeEfficiency / efficiency * 100;
            efficiency = (efficiency * 100) / seeds.size();
            std::cout << "seeds.size():  " << seeds.size() << std::endl;
            std::cout << "layer Zero counter :  " << layerZero << std::endl;
            std::cout << "Avg. efficiency: " << efficiency << "%" << std::endl;
            std::cout << "Avg. fakeEff: " << fakeEfficiency << "%" << std::endl;
            std::cout << "Avg. logFakeEff: " << logFakeEfficiency << "%" << std::endl;
//        }
//    }
    std::string filename = "seed_data.root";
    analyzeSeeds(seeds,filename);
    for (auto s:badSeeds) {
        std::vector<Acts::Vector3> badPos;
        for (auto sl : s.originSourceLinks) {
            badPos.push_back(SA(sl)->
                    localToGlobal(gctx, sl.parameters, Acts::Vector3{0,1,0}));
        }
        for (int b=0;b<badPos.size()-1;b++)
        Acts::GeometryView3D::drawSegment(
                volumeObj,badPos[b],
                badPos[b+1], pConfig);
    }
    volumeObj.write("volumes.obj");

    // Run all configured algorithms and return the appropriate status.
//    return sequencer.run();
      return 0;
} // main
