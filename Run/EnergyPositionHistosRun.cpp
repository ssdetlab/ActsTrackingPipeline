#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXEMagneticField.hpp"

#include <string>
#include <iostream>
/// @brief Run the propagation through 
/// a uniform energy spectrum and record the
/// energy vs position histograms for each layer
int main() {

    Acts::GeometryContext gctx;
    std::string gdmlPath = "lxgeomdump_stave_positron.gdml";
    std::vector<std::string> names{"OPPPSensitive"};

    LUXEGeometry::GeometryOptions gOpt;

    // map (x,y,z) -> (x,y,z)
    auto transformPos = [](const Acts::Vector3& pos) {
        return pos;
    };

    // map (Bx,By,Bz) -> (Bx,By,Bz)
    auto transformBField = [](const Acts::Vector3& field, const Acts::Vector3&) {
        return field;
    };

    const std::vector<unsigned int> bins{5u, 5u, 5u};

    auto BField = LUXEMagneticField::buildLUXEBField(transformPos, transformBField, bins);
    std::cout<<BField.getField(Acts::Vector3{3,1,1}).value()<<std::endl;

    // Build the LUXE detector
    auto positronArmBpr = LUXEGeometry::makeBlueprint(gdmlPath, names, gctx, gOpt);

    auto detector =
        LUXEGeometry::buildLUXEDetector(std::move(positronArmBpr), gctx, gOpt);
    for (auto& vol : detector->rootVolumes()) {
        std::cout<<"Surfaces size: "<<vol->surfaces().size()<<std::endl;
        for (auto& surf : vol->surfaces()) {
            std::cout<<"Surface ID: "<<surf->geometryId()<<std::endl;
            std::cout<<"Surface z transform: "<<surf->center(gctx)[2]<<std::endl;
        }
    }

    return 0;
}