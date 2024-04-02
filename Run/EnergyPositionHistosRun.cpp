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
    std::uniform_real_distribution<> dis1(1.39,2.4);
    std::uniform_real_distribution<> dis2(2.2,5);
    std::uniform_real_distribution<> dis3(5,13);

//    std::gamma_distribution<double> EDis(3, 1.2);
    std::normal_distribution<> pDisP(0.002,0.0018);
    std::normal_distribution<> pDisM(-0.002,0.0018);


    std::vector<LUXENavigator::Measurement> results;
    Acts::ActsScalar m_e = 0.000511;
    std::size_t sourceId = 1;
    for (int i=0;i<200000;i++) {
        Acts::ActsScalar px = (pDisP(gen)+pDisM(gen))/2;
        Acts::ActsScalar pz = (pDisP(gen)+pDisM(gen))/2;
        Acts::ActsScalar E = dis1(gen);
        Acts::ActsScalar py = std::sqrt(std::pow(E,2)-std::pow(m_e,2)-std::pow(px,2)-std::pow(pz,2));
        Acts::ActsScalar p = std::sqrt(std::pow(px,2)+std::pow(py,2)+std::pow(pz,2));
        Acts::ActsScalar theta = std::acos(pz / p);
        Acts::ActsScalar phi = std::atan2(py, px);
        results.push_back(LUXENavigator::createMeasurements(propagator, gctx, magCtx,
                                                            LUXENavigator::makeParameters(p,phi,theta),
                                                            resolutions,sourceId));
        sourceId++;
        if (i%20000==0) std::cout<<"Low E: "<<(i*100)/200000<<"%"<<std::endl;
    };
    saveMeasurementsToFile(results, "low_E_measurements.dat");
    std::vector<LUXENavigator::Measurement> mid_results;
    for (int i=0;i<100000;i++) {
        Acts::ActsScalar px = (pDisP(gen)+pDisM(gen))/2;
        Acts::ActsScalar pz = (pDisP(gen)+pDisM(gen))/2;
        Acts::ActsScalar E = dis2(gen);
        Acts::ActsScalar py = std::sqrt(std::pow(E,2)-std::pow(m_e,2)-std::pow(px,2)-std::pow(pz,2));
        Acts::ActsScalar p = std::sqrt(std::pow(px,2)+std::pow(py,2)+std::pow(pz,2));
        Acts::ActsScalar theta = std::acos(pz / p);
        Acts::ActsScalar phi = std::atan2(py, px);
//        std::cout<<"Initial 4p : "<<px<<" "<<py<<" "<<pz<<std::endl;
//        std::cout<<"Initial dir : "<<phi<<" "<<theta<<std::endl;
        auto mid_res = LUXENavigator::createMeasurements(propagator, gctx, magCtx,
                                                     LUXENavigator::makeParameters(p,phi,theta),
                                                     resolutions,sourceId);
        results.push_back(mid_res);
        mid_results.push_back(mid_res);
        sourceId++;
        if (i%10000==0) std::cout<<"Mid E: "<<(i*100)/100000<<"%"<<std::endl;
    };
    saveMeasurementsToFile(mid_results, "mid_E_measurements.dat");

    std::vector<LUXENavigator::Measurement> high_results;
    for (int i=0;i<20000;i++) {
        Acts::ActsScalar px = (pDisP(gen)+pDisM(gen))/2;
        Acts::ActsScalar pz = (pDisP(gen)+pDisM(gen))/2;
        Acts::ActsScalar E = dis3(gen);
        Acts::ActsScalar py = std::sqrt(std::pow(E,2)-std::pow(m_e,2)-std::pow(px,2)-std::pow(pz,2));
        Acts::ActsScalar p = std::sqrt(std::pow(px,2)+std::pow(py,2)+std::pow(pz,2));
        Acts::ActsScalar theta = std::acos(pz / p);
        Acts::ActsScalar phi = std::atan2(py, px);
//        std::cout<<"Initial 4p : "<<px<<" "<<py<<" "<<pz<<std::endl;
//        std::cout<<"Initial dir : "<<phi<<" "<<theta<<std::endl;
        auto high_res = LUXENavigator::createMeasurements(propagator, gctx, magCtx,
                                                         LUXENavigator::makeParameters(p,phi,theta),
                                                         resolutions,sourceId);
        results.push_back(high_res);
        high_results.push_back(high_res);
        sourceId++;
        if (i%2000==0) std::cout<<"High E: "<<(i*100)/20000<<"%"<<std::endl;
    };
    saveMeasurementsToFile(high_results, "high_E_measurements.dat");
//    for (auto result : results) {
//        std::cout<<result.globalPosition.size()<<std::endl;
//        if (result.globalPosition.size()>1) {
//            for (unsigned int l=0;l<result.globalPosition.size()-1;l++) {
//                Acts::GeometryView3D::drawSegment(
//                        volumeObj,result.globalPosition[l],
//                        result.globalPosition[l+1], pConfig);
//            }
//        }
//    }

    std::string filename = "hist_data_o.root";
    HistogramDatawriter(results,filename,gOpt);

    volumeObj.write("volumes.obj");

    // Run all configured algorithms and return the appropriate status.
//    return sequencer.run();
      return 0;
} // main
