#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXEGeometryIdGenerator.hpp"

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Plugins/Geant4/Geant4SurfaceProvider.hpp"
#include "Acts/Detector/Blueprint.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/detail/BlueprintHelper.hpp"
#include "Acts/Detector/CuboidalContainerBuilder.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"

namespace LUXEGeometry {

std::shared_ptr<Acts::Experimental::LayerStructureBuilder> makeLayerBuilder(std::string gdmlPath,
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
    return LB;
}
std::shared_ptr<const Acts::Experimental::Detector> 
    buildLUXEDetector(std::string gdmlPath, 
        std::vector<std::string> names, 
        Acts::GeometryContext gctx) {
            // Default template parameters are fine
            // when using names as identifiers
            auto layerBuilder1 = makeLayerBuilder(gdmlPath, names, gctx, 3942, 3972);
            auto layerBuilder2 = makeLayerBuilder(gdmlPath, names, gctx, 4042, 4072);
            auto layerBuilder3 = makeLayerBuilder(gdmlPath, names, gctx, 4142, 4172);
            auto layerBuilder4 = makeLayerBuilder(gdmlPath, names, gctx, 4242, 4272);

            // Create the detector
            std::vector<Acts::BinningValue> detectorBins = {Acts::binZ};
            std::vector<Acts::ActsScalar> detectorBounds = {100, 100, 5000};

            // The root node - detector
            auto positronArmBpr = std::make_unique<Acts::Experimental::Blueprint::Node>(
                    "positron_arm", Acts::Transform3::Identity(), Acts::VolumeBounds::eCuboid,
                    detectorBounds, detectorBins);

            std::vector<Acts::ActsScalar> layerBoundaries = {100, 100, 30};

            Acts::Transform3 layer1Transform = Acts::Transform3::Identity() * Acts::Translation3(0., 0., 3962);
            Acts::Transform3 layer2Transform = Acts::Transform3::Identity() * Acts::Translation3(0., 0., 4062);
            Acts::Transform3 layer3Transform = Acts::Transform3::Identity() * Acts::Translation3(0., 0., 4162);
            Acts::Transform3 layer4Transform = Acts::Transform3::Identity() * Acts::Translation3(0., 0., 4262);

            auto layer1Node = std::make_unique<Acts::Experimental::Blueprint::Node>(
                    "layer1", layer1Transform, Acts::VolumeBounds::eCuboid,
                    layerBoundaries, layerBuilder1);
            auto layer2Node = std::make_unique<Acts::Experimental::Blueprint::Node>(
                    "layer2", layer2Transform, Acts::VolumeBounds::eCuboid,
                    layerBoundaries, layerBuilder2);
            auto layer3Node = std::make_unique<Acts::Experimental::Blueprint::Node>(
                    "layer3", layer3Transform, Acts::VolumeBounds::eCuboid,
                    layerBoundaries, layerBuilder3);
            auto layer4Node = std::make_unique<Acts::Experimental::Blueprint::Node>(
                    "layer4", layer4Transform, Acts::VolumeBounds::eCuboid,
                    layerBoundaries, layerBuilder4);

            positronArmBpr->add(std::move(layer1Node));
            positronArmBpr->add(std::move(layer2Node));
            positronArmBpr->add(std::move(layer3Node));
            positronArmBpr->add(std::move(layer4Node));

            positronArmBpr->geoIdGenerator =
                    std::make_shared<Acts::Experimental::LUXEGeometryIdGenerator>(
                            Acts::Experimental::LUXEGeometryIdGenerator::Config{},
                            Acts::getDefaultLogger("RecursiveIdGenerator",
                                                   Acts::Logging::VERBOSE));

            std::cout << "Fill gaps ..." << std::endl;
            // Complete and fill gaps
            Acts::Experimental::detail::BlueprintHelper::fillGaps(*positronArmBpr);
            std::cout << "Filled gaps ..." << std::endl;

            auto detectorBuilder =
                    std::make_shared<Acts::Experimental::CuboidalContainerBuilder>(
                            *positronArmBpr, Acts::Logging::VERBOSE);

            // Detector builder
            Acts::Experimental::DetectorBuilder::Config dCfg;
            dCfg.auxiliary = "*** Test : LUXE detector builder  ***";
            dCfg.name = "LUXE detector from blueprint";
            dCfg.builder = detectorBuilder;
            dCfg.geoIdGenerator = positronArmBpr->geoIdGenerator;

            auto detector = Acts::Experimental::DetectorBuilder(dCfg).construct(gctx);

            return detector;
}

} // namespace LUXEGeometry
