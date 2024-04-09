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

#include "ActsLUXEPipeline/Utils.hpp"

/// @brief Run the propagation through
/// a uniform energy spectrum and record the
/// energy vs position histograms for each layer
int main() {
    using namespace LUXENavigator;

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
    gridOpt.xBins = {-1000,-1, 0.,200, 999,1000.};
    gridOpt.yBins = {1300,1400,1450,1451, 2050.,2649,2650.,2651};
    gridOpt.zBins = {-101,-100, 0.,1, 100, 101.};

    std::string gdmlPath = "lxgeomdump_ip_tracker_positron.gdml";

    std::vector<std::string> staves = {"OPPPSensitive"};
    std::vector<std::string> chamber = {"VCWindowPanel"};
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

    auto positronArmBpr = LUXEGeometry::makeBlueprintLUXE(gdmlPath, staves, gOpt);
    auto detector = LUXEGeometry::buildLUXEDetector(std::move(positronArmBpr), gctx, gOpt);

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
            std::cout<<"Surface x transform: "<<surf->center(gctx)[0]<<std::endl;
            std::cout<<"Surface y transform: "<<surf->center(gctx)[1]<<std::endl;
            std::cout<<"Surface z transform: "<<surf->center(gctx)[2]<<std::endl;
            std::cout<<"Surface bounds: "<<surf->bounds()<<std::endl;
        }
    }

    auto propagator = LUXENavigator::makePropagator<Acts::EigenStepper<>>(detector, BFieldPtr);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> pDisP(0.002,0.0018);
    std::normal_distribution<> pDisM(-0.002,0.0018);
    std::gamma_distribution<double> pzDis(3, 1.2);
    std::uniform_real_distribution<> uni(2.2,2.3);

//    std::vector<Acts::ActsScalar> test_E{0.1,0.3,0.4,0.5,0.6};

    Acts::ActsScalar m_e = 0.000511;
    std::vector<LUXENavigator::Measurement> results;
    std::size_t sourceId = 1;
    int N_events = 100000;
    for (int i=0;i<N_events;i++) {
        Acts::ActsScalar px = (pDisP(gen)+pDisM(gen))/2;
        Acts::ActsScalar pz = (pDisP(gen)+pDisM(gen))/2;
        Acts::ActsScalar py = pzDis(gen)+1;
//        Acts::ActsScalar py = uni(gen);
        Acts::ActsScalar p = std::sqrt(std::pow(px,2)+std::pow(py,2)+std::pow(pz,2));
        Acts::ActsScalar E = std::hypot(p,m_e);
        Acts::ActsScalar theta = std::acos(pz / p);
        Acts::ActsScalar phi = std::atan2(py, px);
        results.push_back(LUXENavigator::createMeasurements(propagator, gctx, mctx,
                                                            LUXENavigator::makeParameters(p,phi,theta),sourceId));
        sourceId++;
        if (i%(N_events/10)==0) {
            std::cout<<"Completed: "<<(i*100)/N_events<<"%"<<std::endl;
        }
    };

    saveMeasurementsToFile(results, "100k_measurements.dat");



    return 0;
} // main
