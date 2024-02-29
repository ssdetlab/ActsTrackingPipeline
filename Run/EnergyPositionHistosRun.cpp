#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXEMagneticField.hpp"

#include <string>
#include <iostream>
/// @brief Run the propagation through 
/// a uniform energy spectrum and record the
/// energy vs position histograms for each layer
int main() {
    // Get the mockup Geant4 world
    auto [world, names, gdmlPath] = LUXEGeometry::MockupGeant4World();

    Acts::GeometryContext gctx;
    std::string gdmlPath2 = "lxgeomdump_staves.gdml";
    std::vector<std::string> names2{"OPPPSensitive"};

    // Build the LUXE detector
    auto detector =
        LUXEGeometry::buildLUXEDetector(gdmlPath2, names2, gctx);


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
        // linear in r and z so interpolation should be exact
        return Acts::Vector3(3 * x * z / r5, 3 * y * z / r5,
                             (3 * std::pow(z, 2) - std::pow(r, 2)) / r5);
    };

    const std::vector<unsigned int> bins{5u, 5u, 5u};

    auto BField = LUXEMagneticField::buildLUXEBField(transformPos, transformBField, bFieldValue, bins);
    std::cout<<BField.getField(Acts::Vector3{3,1,1}).value()<<std::endl;

    return 0;
}