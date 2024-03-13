#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXEGeometryIdGenerator.hpp"

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Plugins/Geant4/Geant4SurfaceProvider.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/detail/BlueprintHelper.hpp"
#include "Acts/Detector/CuboidalContainerBuilder.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"
#include "Acts/Definitions/Units.hpp"

#include <vector>

namespace LUXEGeometry {
using namespace Acts::UnitLiterals;
using layerBuilderVector = std::vector<std::shared_ptr<Acts::Experimental::LayerStructureBuilder>>;
using layerNodeVector = std::vector<std::unique_ptr<Acts::Experimental::Blueprint::Node>>;

std::shared_ptr<Acts::Experimental::LayerStructureBuilder> makeLayerBuilder(std::string gdmlPath,
                                                                            std::vector<std::string> names,
                                                                            Acts::GeometryContext gctx,
                                                                            Acts::ActsScalar zMin,
                                                                            Acts::ActsScalar zMax) {
    auto spCfg = Acts::Experimental::Geant4SurfaceProvider<1>::Config();
    spCfg.gdmlPath = gdmlPath;
    spCfg.surfacePreselector =
            std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(names,
                                                                                true);

    auto kdtDOpt = Acts::Experimental::Geant4SurfaceProvider<1>::kdtOptions();
    kdtDOpt.range[0].set(zMin, zMax);
    kdtDOpt.binningValues = {Acts::BinningValue::binZ};

    auto SP = std::make_shared<Acts::Experimental::Geant4SurfaceProvider<1>>(
            spCfg, kdtDOpt, false);

    auto lbCfg = Acts::Experimental::LayerStructureBuilder::Config();
    lbCfg.surfacesProvider = SP;
    auto LB =
            std::make_shared<Acts::Experimental::LayerStructureBuilder>(lbCfg);
    return LB;
}

blueprintPtr makeBlueprint(std::string gdmlPath,
              std::vector<std::string> names,
              Acts::GeometryContext gctx,
              LUXEGeometry::GeometryOptions gOpt) {
        size_t numLayers = gOpt.layerZPositions.size();

        layerBuilderVector layerBuilders(numLayers);

        std::vector<Acts::BinningValue> detectorBins = {Acts::binZ};

        auto positronArmBpr = std::make_unique<Acts::Experimental::Blueprint::Node>(
                "positron_arm", Acts::Transform3::Identity(), Acts::VolumeBounds::eCuboid,
                gOpt.detectorBounds, detectorBins);

        std::vector<Acts::Transform3> layerTransforms(numLayers);
        layerNodeVector layerNodes(numLayers);

        for (int i = 0; i<static_cast<int>(numLayers);i++) {
            std::cout<<"BOUNDS: "<< gOpt.layerZPositions[i] - gOpt.deltaZ- 1_mm<<" "
                                <<gOpt.layerZPositions[i] + 1_mm<<std::endl;

            layerBuilders[i] = makeLayerBuilder(gdmlPath, names, gctx,
                                                gOpt.layerZPositions[i] - gOpt.deltaZ- 1_mm,
                                                gOpt.layerZPositions[i] + 1_mm);
            layerTransforms[i] = Acts::Transform3::Identity() *
                                 Acts::Translation3(0., 0., gOpt.layerZPositions[i]);


            layerNodes[i] = std::make_unique<Acts::Experimental::Blueprint::Node>(
                    "layer"+std::to_string(i), layerTransforms[i], Acts::VolumeBounds::eCuboid,
                    gOpt.layerBounds, layerBuilders[i]);

            positronArmBpr->add(std::move(layerNodes[i]));
        }
        return positronArmBpr;
}


std::shared_ptr<const Acts::Experimental::Detector>
    buildLUXEDetector(const blueprintPtr detectorBpr,
        Acts::GeometryContext gctx,
        LUXEGeometry::GeometryOptions gOpt) {

            detectorBpr->geoIdGenerator =
                    std::make_shared<Acts::Experimental::LUXEGeometryIdGenerator>(
                            Acts::Experimental::LUXEGeometryIdGenerator::Config{},
                            Acts::getDefaultLogger("RecursiveIdGenerator",
                                                   Acts::Logging::VERBOSE));

            std::cout << "Fill gaps ..." << std::endl;
            // Complete and fill gaps
            Acts::Experimental::detail::BlueprintHelper::fillGaps(*detectorBpr);
            std::cout << "Filled gaps ..." << std::endl;

            auto detectorBuilder =
                    std::make_shared<Acts::Experimental::CuboidalContainerBuilder>(
                            *detectorBpr, Acts::Logging::VERBOSE);

            // Detector builder
            Acts::Experimental::DetectorBuilder::Config dCfg;
            dCfg.auxiliary = "*** Test : LUXE detector builder  ***";
            dCfg.name = "LUXE detector from blueprint";
            dCfg.builder = detectorBuilder;
            dCfg.geoIdGenerator = detectorBpr->geoIdGenerator;

            auto detector = Acts::Experimental::DetectorBuilder(dCfg).construct(gctx);

            return detector;
}

} // namespace LUXEGeometry
