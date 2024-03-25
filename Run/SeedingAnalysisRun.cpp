#include "ActsLUXEPipeline/Sequencer.hpp"

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

#include "ActsLUXEPipeline/Utils.hpp"

/// @brief Run the propagation through
/// a uniform energy spectrum and record the
/// energy vs position histograms for each layer
int main() {
    using namespace LUXENavigator;
    Acts::Logging::Level logLevel = Acts::Logging::VERBOSE;

    // setup the sequencer first w/ config derived from options
    Sequencer::Config seqCfg;
    seqCfg.events = 10;
    seqCfg.numThreads = -1;
    Sequencer sequencer(seqCfg);

//    LUXEROOTReader::LUXEROOTSimDataReader::Config readerCfg
//        = LUXEROOTReader::defaultSimConfig();
//    readerCfg.dataCollection = "SourceLink";
//    std::string pathToDir = "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_dataInRootFormat/SignalNextTrial_e1gpc_10.0_1";
    // map (x,y,z) -> (x,y,z)

    auto transformPos = [](const Acts::Vector3& pos) {
        LUXEGeometry::GeometryOptions gOpt;
        for (int i=0;i<3;i++) {
            if (pos[i]<=gOpt.MagneticFieldBounds[i].first ||
                pos[i]>gOpt.MagneticFieldBounds[i].second) {
                return Acts::Vector3{0,1300,0};
            }
        }
        return pos;
    };

    // map (Bx,By,Bz) -> (Bx,By,Bz)
    auto transformBField = [](const Acts::Vector3& field, const Acts::Vector3&) {
        return field;
    };

    LUXEMagneticField::GridOptions gridOpt;
    gridOpt.bins = {14u, 1400u, 14u};
    gridOpt.limits = {std::make_pair(-600,800),
                      std::make_pair(1250,2850),
                      std::make_pair(-600,800)};

    auto BField = LUXEMagneticField::buildLUXEBField(transformPos, transformBField, gridOpt);
    auto BFieldPtr = std::make_shared<LUXEMagneticField::BField_t>(BField);

    // Build the LUXE detector
    std::string gdmlPath = "lxgeomdump_ip_tracker_positron.gdml";

    std::vector<std::string> staves = {"OPPPSensitive"};
    std::vector<std::string> chamber = {"VCWindowPanel"};
    Acts::GeometryContext gctx;
    Acts::MagneticFieldContext magCtx;
    LUXEGeometry::GeometryOptions gOpt;
    auto magneticChamberBpr = LUXEGeometry::makeBlueprintMagneticChamber(gdmlPath, chamber, gOpt);
    auto positronArmBpr = LUXEGeometry::makeBlueprintPositron(gdmlPath, staves, gOpt);
    positronArmBpr->add(std::move(magneticChamberBpr));
    auto detector = LUXEGeometry::buildLUXEDetector(std::move(positronArmBpr), gctx, gOpt);

    MeasurementResolution resPixel = {MeasurementType::eLoc01,
                                      {gOpt.chipSizeX,
                                       gOpt.chipSizeY}};
    std::vector<std::pair<Acts::GeometryIdentifier,MeasurementResolution>> m;
    Acts::ViewConfig pConfig = Acts::s_viewSensitive;
    Acts::ObjVisualization3D volumeObj;
    for (auto& vol : detector->rootVolumes()) {
        std::cout<<"Surfaces size: "<<vol->surfaces().size()<<std::endl;
        std::cout<<"Volume Bounds: "<<vol->volumeBounds()<<std::endl;
        std::cout<<"Volume Transformation: "<<vol->transform().translation()<<std::endl;
//        Acts::GeometryView3D::drawDetectorVolume(
//                volumeObj, *(vol), gctx,
//                Acts::Transform3::Identity(), pConfig);
        for (auto& surf : vol->surfaces()) {
            std::cout<<"Assigning resolution to surface ID: "<<surf->geometryId()<<std::endl;
//            if (vol->geometryId().volume()!=1) {
                Acts::GeometryView3D::drawSurface(
                        volumeObj, *(surf), gctx,
                        Acts::Transform3::Identity(), pConfig);
//            }
            m.push_back(std::make_pair(surf->geometryId(),resPixel));
            std::cout<<"Surface x transform: "<<surf->center(gctx)[0]<<std::endl;
            std::cout<<"Surface y transform: "<<surf->center(gctx)[1]<<std::endl;
            std::cout<<"Surface z transform: "<<surf->center(gctx)[2]<<std::endl;
            std::cout<<"Surface bounds: "<<surf->bounds()<<std::endl;
        }
    }
    MeasurementResolutionMap resolutions = m;

    auto propagator = LUXENavigator::makePropagator<Acts::EigenStepper<>>(detector, BFieldPtr);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> pDisP(0.002,0.0018);
    std::normal_distribution<> pDisM(-0.002,0.0018);
    std::gamma_distribution<double> pzDis(3, 1.2);
    std::uniform_real_distribution<> uni(8,12);
// lower limit 0.32 GeV upper limit 3.4 GeV

//    std::vector<Acts::ActsScalar> test_E{0.1,0.3,0.4,0.5,0.6};

    Acts::ActsScalar m_e = 0.000511;
    std::vector<LUXENavigator::Measurements> results;
    std::size_t sourceId = 1;

    for (int i=0;i<5000;i++) {
        Acts::ActsScalar px = (pDisP(gen)+pDisM(gen))/2;
        Acts::ActsScalar pz = (pDisP(gen)+pDisM(gen))/2;
        Acts::ActsScalar py = pzDis(gen)+1;
//        Acts::ActsScalar py = uni(gen);
        Acts::ActsScalar p = std::sqrt(std::pow(px,2)+std::pow(py,2)+std::pow(pz,2));
        Acts::ActsScalar E = std::hypot(p,m_e);
        std::cout<<"Initial 4p : "<<px<<" "<<py<<" "<<pz<<" "<<E<<std::endl;
        Acts::ActsScalar theta = std::acos(pz / p);
        Acts::ActsScalar phi = std::atan2(py, px);
        std::cout<<"Initial dir : "<<phi<<" "<<theta<<std::endl;
        results.push_back(LUXENavigator::createMeasurements(propagator, gctx, magCtx,
                                                            LUXENavigator::makeParameters(p,phi,theta),
                                                            resolutions,sourceId));
        sourceId++;
    };

    SimpleSourceLink::SurfaceAccessor SA{*detector};
    std::vector<SimpleSourceLink> sl4Seeding;
    for (auto& result:results) {
        for (unsigned int l=0;l<result.truthParameters.size();l++) {
            sl4Seeding.push_back(result.sourceLinks[l]);
        } //static_cast<float>
//        for (unsigned int l=1;l<result.fullTrack.size()-3;l++) {
//            if (result.fullTrack[l][1]<6000) {
//                Acts::GeometryView3D::drawSegment(
//                        volumeObj,result.fullTrack[l],
//                        result.fullTrack[l+1], pConfig);
//            }
//        }
    }
    std::vector<LUXETrackFinding::Seed> trueSeeds;
    std::vector<SimpleSourceLink> seedSLs;

    for (size_t s=0 ; s < sl4Seeding.size(); s++) {
        seedSLs = {sl4Seeding[s]};
        size_t i=0;
        if (s+1 == sl4Seeding.size()) {
            LUXETrackFinding::Seed seed{sl4Seeding[s],0,{0},seedSLs,Acts::Vector4{0,0,0,0}};
            trueSeeds.push_back(seed);
            continue;
        }
        while (sl4Seeding[s].eventId==sl4Seeding[s+1+i].eventId && s+i < sl4Seeding.size()-1) {
            seedSLs.push_back(sl4Seeding[s+1+i]);
            i++;
        }
        s+=i;
        LUXETrackFinding::Seed seed{sl4Seeding[s],0,{0},seedSLs,Acts::Vector4{0,0,0,0}};
        trueSeeds.push_back(seed);
        if (sl4Seeding[s].eventId+1<sl4Seeding[s+1].eventId) {
            LUXETrackFinding::Seed seed{sl4Seeding[s],0,{0},{},Acts::Vector4{0,0,0,0}};
            trueSeeds.push_back(seed);
        }
    }

    std::vector<LUXETrackFinding::Seed> seeds =
            LUXETrackFinding::LUXEPathSeeder(gctx, gOpt, detector, sl4Seeding,
                                             "/Users/alonlevi/CLionProjects/LUXEPipeline/build");

    std::cout<<"Seeds Produced : "<<seeds.size()<<std::endl;

    SimpleSourceLink::SurfaceAccessor SAseed{*detector};
    int count = 1;
    pConfig.color = {0,250,0};
    for (auto& seed : seeds) {
        std::cout<<"source links in seed "<<count<<": "<<seed.sourceLinks.size()<<std::endl;
        count++;
        std::vector<Acts::Vector3> seedPos;
        for (auto& sl : seed.sourceLinks) {
            seedPos.push_back(SAseed(sl)->
                    localToGlobal(gctx, sl.parameters, Acts::Vector3{0,1,0}));
        }
    }
//    for (auto sl:sl4Seeding) {
//        std::cout<<"sl4seeding event "<<sl.eventId<<std::endl;
//    }
//    int trueSeedCount = 0;
//    for (auto seed : trueSeeds) {
//        trueSeedCount++;
//        std::cout<<"True Seed #: "<<trueSeedCount<<std::endl;
//        for (auto sl : seed.sourceLinks) {
//            std::cout<<sl.eventId<<std::endl;
//        }
//    }
//    int seedCount = 0;
//    for (auto seed : seeds) {
//        seedCount++;
//        std::cout<<"Seed #: "<<seedCount<<std::endl;
//        for (auto sl : seed.sourceLinks) {
//            std::cout<<sl.eventId<<std::endl;
//        }
//    }


    std::string filename = "seed_data.root";
    analyzeSeeds(seeds,filename);
    volumeObj.write("volumes.obj");

    // Run all configured algorithms and return the appropriate status.
//    return sequencer.run();
      return 0;
} // main
