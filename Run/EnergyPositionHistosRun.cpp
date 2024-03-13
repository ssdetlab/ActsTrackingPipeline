#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXEGeometryConstraints.hpp"
#include "ActsLUXEPipeline/LUXEMagneticField.hpp"
#include "ActsLUXEPipeline/LUXENavigator.hpp"
#include "ActsLUXEPipeline/LUXEMeasurementsCreator.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include <string>
#include <iostream>
/// @brief Run the propagation through 
/// a uniform energy spectrum and record the
/// energy vs position histograms for each layer
int main() {
    using namespace LUXENavigator;
    Acts::GeometryContext gctx;
    std::string gdmlPath = "lxgeomdump_stave_positron.gdml";
    std::vector<std::string> names{"OPPPSensitive"};
    LUXEGeometry::GeometryOptions gOpt;

    Acts::MagneticFieldContext magCtx;
    // map (x,y,z) -> (x,y,z)
    auto transformPos = [](const Acts::Vector3& pos) {
        return pos;
    };

    // map (Bx,By,Bz) -> (Bx,By,Bz)
    auto transformBField = [](const Acts::Vector3& field, const Acts::Vector3&) {
        return field;
    };

    const std::vector<unsigned int> bins{10u, 10u, 10u};

    auto BField = LUXEMagneticField::buildLUXEBField(transformPos, transformBField, bins);
    std::cout<<BField.getField(Acts::Vector3{3,1,1}).value()<<std::endl;
    auto BFieldPtr = std::make_shared<LUXEMagneticField::BField_t>(BField);

    // Build the LUXE detector
    auto positronArmBpr = LUXEGeometry::makeBlueprint(gdmlPath, names, gctx, gOpt);

    auto testParams = LUXENavigator::makeParameters();

    MeasurementResolution resPixel = {MeasurementType::eLoc01,
                                      {LUXEGeometry::chipSizeX,
                                       LUXEGeometry::chipSizeY}};

    Acts::ViewConfig pConfig = Acts::s_viewSensitive;

    Acts::ObjVisualization3D volumeObj;
    std::vector<std::pair<Acts::GeometryIdentifier,MeasurementResolution>> m;
    auto detector =
            LUXEGeometry::buildLUXEDetector(std::move(positronArmBpr), gctx, gOpt);

    for (auto& vol : detector->rootVolumes()) {
        std::cout<<"Surfaces size: "<<vol->surfaces().size()<<std::endl;
        Acts::GeometryView3D::drawDetectorVolume(
                volumeObj, *(vol), gctx,
                Acts::Transform3::Identity(), pConfig);
        for (auto& surf : vol->surfaces()) {
            std::cout<<"Assigning resolution to surface ID: "<<surf->geometryId()<<std::endl;
            Acts::GeometryView3D::drawSurface(
                    volumeObj, *(surf), gctx,
                    Acts::Transform3::Identity(), pConfig);
            m.push_back(std::make_pair(surf->geometryId(),resPixel));
            std::cout<<"Surface x transform: "<<surf->center(gctx)[0]<<std::endl;
            std::cout<<"Surface y transform: "<<surf->center(gctx)[1]<<std::endl;
            std::cout<<"Surface bounds: "<<surf->bounds()<<std::endl;
        }
    }



    MeasurementResolutionMap resolutions = m;

    auto propagator = LUXENavigator::makePropagator<Acts::EigenStepper<>>(detector, BFieldPtr);

    auto test = LUXENavigator::createMeasurements(propagator, gctx, magCtx, testParams, resolutions);

//    Acts::GeometryView3D::drawSegment(
//            volumeObj, *(surf), gctx,
//            Acts::Transform3::Identity(), pConfig);
    for (auto& sl:test.sourceLinks) {
        std::cout<<sl.parameters<<std::endl;
    }
    volumeObj.write("volumes.obj");

    return 0;
}