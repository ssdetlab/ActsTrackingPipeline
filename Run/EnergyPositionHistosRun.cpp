#include "ActsLUXEPipeline/LUXEGeometry.hpp"

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

    // Build the LUXE detector
    auto positronArmBpr = LUXEGeometry::makeBlueprint(gdmlPath, names, gctx, gOpt);

    auto detector =
        LUXEGeometry::buildLUXEDetector(positronArmBpr, gctx, gOpt);
    for (auto& vol : detector->rootVolumes()) {
        std::cout<<"Surfaces size: "<<vol->surfaces().size()<<std::endl;
        for (auto& surf : vol->surfaces()) {
            std::cout<<"Surface ID: "<<surf->geometryId()<<std::endl;
            std::cout<<"Surface z transform: "<<surf->center(gctx)[2]<<std::endl;
        }
    }

    return 0;
}