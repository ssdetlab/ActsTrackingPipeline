#include "ActsLUXEPipeline/LUXEGeometry.hpp"

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Plugins/Geant4/Geant4SurfaceProvider.hpp"

namespace LUXEGeometry {

std::shared_ptr<const Acts::Experimental::Detector> 
    buildLUXEDetector(std::string gdmlPath, 
        std::vector<std::string> names, 
        Acts::GeometryContext gctx) {
            // Default template parameters are fine
            // when using names as identifiers
            auto spFullCfg = Acts::Experimental::Geant4SurfaceProvider<>::Config();
            spFullCfg.gdmlPath = gdmlPath;
            spFullCfg.surfacePreselector =
                std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(names,
                    true);
    
            auto spFull = std::make_shared<Acts::Experimental::Geant4SurfaceProvider<>>(
                spFullCfg, Acts::Experimental::Geant4SurfaceProvider<>::kdtOptions(),
                false);
    
            auto lbFullCfg = Acts::Experimental::LayerStructureBuilder::Config();
            lbFullCfg.surfacesProvider = spFull;
    
            auto lbFull =
                std::make_shared<Acts::Experimental::LayerStructureBuilder>(lbFullCfg);
    
            auto [sFull, vFull, suFull, vuFull] = lbFull->construct(gctx);
    
            assert(sFull.size() == names.size());

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
