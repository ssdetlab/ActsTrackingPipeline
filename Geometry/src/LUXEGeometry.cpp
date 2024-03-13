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

std::unique_ptr<Acts::Experimental::Blueprint::Node>
makeBlueprintPositron(
    const std::string& gdmlPath,
    const std::vector<std::string>& names,
    const LUXEGeometry::GeometryOptions& gOpt) {
        // Read the gdml file and get the world volume
        G4GDMLParser parser;
        parser.Read(gdmlPath, false);
        auto world = parser.GetWorldVolume();

        // Transformation has to be applied to the world
        // because the z-spacing of the detector
        // and the Acts track parametrization
        // do not work together
        G4Transform3D toWorld(gOpt.g4WorldRotation,
            G4ThreeVector(0, 0, 0));

        std::size_t numLayers = gOpt.layerZPositions.size();
        std::vector<Acts::BinningValue> detectorBins = {Acts::binZ};

        // We adopt the convention of first
        // rotating and then translating
        // (order matters here)
        Acts::Transform3 positronArmTransform = Acts::Transform3::Identity();
        positronArmTransform.rotate(gOpt.actsWorldRotation);
        positronArmTransform.translate(gOpt.postironArmTranslation);

        // Create the root node of the blueprint
        auto positronArmBpr =
            std::make_unique<Acts::Experimental::Blueprint::Node>(
            "positronArm", positronArmTransform,
            Acts::VolumeBounds::eCuboid,
            gOpt.postironArmBounds, detectorBins);

        // Iterate over the layers and create
        // the child nodes
        for (std::size_t i = 0; i < numLayers; i++) {
            auto zBounds = std::make_tuple(
                gOpt.layerZPositions[i] - gOpt.deltaZ,
                gOpt.layerZPositions[i] + gOpt.deltaZ);

            // As the volumes are already rotated,
            // the selection has to happen along the y-axis
            auto layerBuilder = makeLayerBuilder(
                world, toWorld, names, {zBounds}, {Acts::binY});

            // Convention is that the transformations
            // are with respect to the global frame
            Acts::Vector3 layerTranslation =
                Acts::Vector3(
                    gOpt.postironArmTranslation.x(),
                    gOpt.postironArmTranslation.y(),
                    gOpt.layerZPositions[i]);

            Acts::Transform3 layerTransform = Acts::Transform3::Identity();
            layerTransform.rotate(gOpt.actsWorldRotation);
            layerTransform.translate(layerTranslation);

            auto layerNode = std::make_unique<Acts::Experimental::Blueprint::Node>(
                "layer" + std::to_string(i), layerTransform,
                Acts::VolumeBounds::eCuboid, gOpt.layerBounds,
                layerBuilder);

            positronArmBpr->add(std::move(layerNode));
        }

        return positronArmBpr;
}


std::shared_ptr<const Acts::Experimental::Detector>
    buildLUXEDetector(
        const std::unique_ptr<Acts::Experimental::Blueprint::Node>
            detectorBpr,
        const Acts::GeometryContext& gctx,
        const LUXEGeometry::GeometryOptions& gOpt) {
            detectorBpr->geoIdGenerator =
                std::make_shared<LUXEGeometryIdGenerator>(
                    LUXEGeometryIdGenerator::Config{},
                    Acts::getDefaultLogger("GeoIdGenerator",
                        Acts::Logging::VERBOSE));

        // Complete and fill gaps
        Acts::Experimental::detail::BlueprintHelper::fillGaps(*detectorBpr, false);

        auto detectorBuilder =
            std::make_shared<Acts::Experimental::CuboidalContainerBuilder>(
                *detectorBpr, Acts::Logging::VERBOSE);

        // Detector builder
        Acts::Experimental::DetectorBuilder::Config dCfg;
        dCfg.auxiliary = "LUXE detector builder";
        dCfg.name = "LUXE detector from blueprint";
        dCfg.builder = detectorBuilder;
        dCfg.geoIdGenerator = detectorBpr->geoIdGenerator;

        auto detector = Acts::Experimental::DetectorBuilder(dCfg).construct(gctx);

        return detector;
}

} // namespace LUXEGeometry
