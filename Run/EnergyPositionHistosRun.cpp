#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXEMagneticField.hpp"

#include "Acts/Definitions/Algebra.hpp"
/// @brief Run the propagation through 
/// a uniform energy spectrum and record the
/// energy vs position histograms for each layer
int main() {
    // Get the mockup Geant4 world
    auto [world, names, gdmlPath] = LUXEGeometry::MockupGeant4World();

    Acts::GeometryContext gctx;

    // Build the LUXE detector
    auto detector =
        LUXEGeometry::buildLUXEDetector(gdmlPath, names, gctx);

    auto BField = LUXEMagneticField::buildLUXEBField();
    std::cout<<BField.getField(Acts::Vector3{3,1,1}).value()<<std::endl;

    return 0;
}