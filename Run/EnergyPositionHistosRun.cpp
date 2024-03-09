#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXEMagneticField.hpp"

#include <string>
#include <iostream>
/// @brief Run the propagation through 
/// a uniform energy spectrum and record the
/// energy vs position histograms for each layer
int main() {
    // Get the mockup Geant4 world

    Acts::GeometryContext gctx;
    std::string gdmlPath = "lxgeomdump_stave_positron.gdml";
    std::vector<std::string> names{"OPPPSensitive"};

    // Build the LUXE detector
    auto detector =
        LUXEGeometry::buildLUXEDetector(gdmlPath, names, gctx);


    // map (x,y,z) -> (x,y,z)
    auto transformPos = [](const Acts::Vector3& pos) {
        return pos;
    };

    // map (Bx,By,Bz) -> (Bx,By,Bz)
    auto transformBField = [](const Acts::Vector3& field, const Acts::Vector3&) {
        return field;
    };

    auto bFieldValue = [](const std::array<double, 3> &v) {
        double x = v.at(0);
        double y = v.at(1);
        double z = v.at(2);
        double r = std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2));
        double r5 = std::pow(r,5);

        return Acts::Vector3(3 * x * z / r5, 3 * y * z / r5,
                             (3 * std::pow(z, 2) - std::pow(r, 2)) / r5);
    };

    const std::vector<unsigned int> bins{5u, 5u, 5u};

    auto BField = LUXEMagneticField::buildLUXEBField(transformPos, transformBField, bFieldValue, bins);

    for (auto& vol : detector->rootVolumes()) {
        std::cout<<"Surfaces size: "<<vol->surfaces().size()<<std::endl;
        for (auto& surf : vol->surfaces()) {
            std::cout<<"Surface ID: "<<surf->geometryId()<<std::endl;
            std::cout<<"Surface z transform: "<<surf->center(gctx)[2]<<std::endl;
        }
    }

    return 0;
}