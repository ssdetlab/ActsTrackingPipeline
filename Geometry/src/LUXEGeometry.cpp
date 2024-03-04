#include "ActsLUXEPipeline/LUXEGeometry.hpp"

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Plugins/Geant4/Geant4SurfaceProvider.hpp"
namespace LUXEGeometry {

Acts::Experimental::InternalStructure selectSurfaces(std::string gdmlPath,
                                                     std::vector<std::string> names,
                                                     Acts::GeometryContext gctx,
                                                     Acts::ActsScalar zMin,
                                                     Acts::ActsScalar zMax) {
    auto spCfg = Acts::Experimental::Geant4SurfaceProvider<3>::Config();
    spCfg.gdmlPath = gdmlPath;
    spCfg.surfacePreselector =
            std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(names,
                                                                                true);

    auto kdtDOpt = Acts::Experimental::Geant4SurfaceProvider<3>::kdtOptions();
    kdtDOpt.range = Acts::RangeXD<3, Acts::ActsScalar>();
    kdtDOpt.range[0].set(zMin, zMax);
    kdtDOpt.binningValues = {Acts::BinningValue::binZ};

    auto SP = std::make_shared<Acts::Experimental::Geant4SurfaceProvider<3>>(
            spCfg, kdtDOpt, false);

    auto lbCfg = Acts::Experimental::LayerStructureBuilder::Config();
    lbCfg.surfacesProvider = SP;
    auto LB =
            std::make_shared<Acts::Experimental::LayerStructureBuilder>(lbCfg);

    return LB->construct(gctx);
}
std::shared_ptr<const Acts::Experimental::Detector> 
    buildLUXEDetector(std::string gdmlPath, 
        std::vector<std::string> names, 
        Acts::GeometryContext gctx) {
            // Default template parameters are fine
            // when using names as identifiers
            auto [s1, v1, su1, vu1] = selectSurfaces(gdmlPath, names, gctx, 3940, 3972);
            auto [s2, v2, su2, vu2] = selectSurfaces(gdmlPath, names, gctx, 4040, 4072);
            auto [s3, v3, su3, vu3] = selectSurfaces(gdmlPath, names, gctx, 4140, 4172);
            auto [s4, v4, su4, vu4] = selectSurfaces(gdmlPath, names, gctx, 4240, 4272);
            std::cout<<"s1 Size: "<<s1.size()<<std::endl;
            std::cout<<"s2 Size: "<<s2.size()<<std::endl;
            std::cout<<"s3 Size: "<<s3.size()<<std::endl;
            std::cout<<"s4 Size: "<<s4.size()<<std::endl;

            // Create the detector
            auto bounds = std::make_unique<Acts::CuboidVolumeBounds>(3, 3, 3);
            auto volume = Acts::Experimental::DetectorVolumeFactory::construct(
                Acts::Experimental::defaultPortalAndSubPortalGenerator(), gctx,
                "volume", Acts::Transform3::Identity(), std::move(bounds), {}, {},
                Acts::Experimental::tryNoVolumes(),
                Acts::Experimental::tryAllPortalsAndSurfaces());
            volume->assignGeometryId(1);
            auto detector = Acts::Experimental::Detector::makeShared(
                "detector", {volume}, Acts::Experimental::tryRootVolumes());
            
            return detector;
}

} // namespace LUXEGeometry
