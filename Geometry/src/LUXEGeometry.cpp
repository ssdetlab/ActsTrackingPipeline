#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXEGeometryIdGenerator.hpp"

#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/detail/BlueprintHelper.hpp"
#include "Acts/Detector/CuboidalContainerBuilder.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"

#include <vector>

namespace LUXEGeometry {

std::unique_ptr<Acts::Experimental::Blueprint::Node> 
makeBlueprintLUXE(
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

        std::vector<Acts::BinningValue> worldBins = {Acts::binZ};
        std::vector<Acts::BinningValue> trackerBins = {Acts::binX};

        Acts::Transform3 trackerTransform = Acts::Transform3::Identity();
        trackerTransform.rotate(gOpt.actsWorldRotation);
        trackerTransform.translate(gOpt.trackerTranslation);

        // Create the tracker node of the blueprint
        auto trackerBP = 
            std::make_unique<Acts::Experimental::Blueprint::Node>(
            "Tracker", trackerTransform,
            Acts::VolumeBounds::eCuboid,
            gOpt.trackerBounds, trackerBins);

        // We adopt the convention of first 
        // rotating and then translating
        // (order matters here)
        Acts::Transform3 positronArmTransform = Acts::Transform3::Identity();
        positronArmTransform.prerotate(gOpt.actsWorldRotation);
        positronArmTransform.translate(gOpt.postironArmTranslation);

        Acts::Transform3 electronArmTransform = Acts::Transform3::Identity();
        electronArmTransform.rotate(gOpt.actsWorldRotation);
        electronArmTransform.translate(gOpt.electronArmTranslation);

        // Create arm nodes of the blueprint
        auto positronArmBP = 
            makeBlueprintArm(
                world, toWorld, "Positron", positronArmTransform, 
                gOpt.armBounds, names, 
                gOpt.layerZPositions, gOpt.layerBounds);

        auto electronArmBP =
            makeBlueprintArm(
                world, toWorld, "Electron", electronArmTransform, 
                gOpt.armBounds, names, 
                gOpt.layerZPositions, gOpt.layerBounds);

        trackerBP->add(std::move(positronArmBP));
        // trackerBP->add(std::move(electronArmBP));

        return trackerBP;
};

std::unique_ptr<Acts::Experimental::Blueprint::Node> 
makeBlueprintArm(
    const G4VPhysicalVolume* world,
    const G4Transform3D& toWorld,
    const std::string& armName,
    const Acts::Transform3& armTransform,
    const std::vector<Acts::ActsScalar>& armBounds,
    const std::vector<std::string>& chipNames,
    const std::vector<Acts::ActsScalar>& layerZPositions,
    const std::vector<Acts::ActsScalar>& layerBounds) {
        std::size_t numLayers = layerZPositions.size();
        std::vector<Acts::BinningValue> armBins = {Acts::binZ};

        // Create arm nodes of the blueprint
        auto armBP = 
            std::make_unique<Acts::Experimental::Blueprint::Node>(
            armName + "Arm", armTransform,
            Acts::VolumeBounds::eCuboid,
            armBounds, armBins);

        // Iterate over the layers and create 
        // the child nodes
        for (std::size_t i = 0; i < numLayers; i++) {
            // Layer bounds
            auto zBounds = std::make_tuple(
                layerZPositions[i] - layerBounds[2],
                layerZPositions[i] + layerBounds[2]);
            auto xBounds = std::make_tuple(
                armTransform.translation().x() - armBounds[0], 
                armTransform.translation().x() + armBounds[0]);

            // As the volumes are already rotated, 
            // the selection has to happen along the y-axis
            auto layerBuilder = makeLayerBuilder<2>(
                world, toWorld, chipNames, 
                {xBounds, zBounds}, 
                {Acts::binX, Acts::binY});

            // Convention is that the transformations
            // are with respect to the global frame
            Acts::Vector3 layerTranslation =
                Acts::Vector3(
                    armTransform.translation().x(), 
                    layerZPositions[i],
                    armTransform.translation().z()); 

            Acts::Transform3 layerTransform = Acts::Transform3::Identity();
            layerTransform.rotate(armTransform.rotation());
            layerTransform.pretranslate(layerTranslation);

            auto layerNode = 
                std::make_unique<Acts::Experimental::Blueprint::Node>(
                    "layer" + armName + std::to_string(i), layerTransform, 
                    Acts::VolumeBounds::eCuboid, layerBounds, 
                    layerBuilder);
    
            armBP->add(std::move(layerNode));
        }

        return armBP;
};

std::shared_ptr<const Acts::Experimental::Detector>
    buildLUXEDetector(
        const std::unique_ptr<Acts::Experimental::Blueprint::Node> 
            detectorBpr,
        const Acts::GeometryContext& gctx,
        const LUXEGeometry::GeometryOptions& gOpt) {
            // detectorBpr->geoIdGenerator =

        // Complete and fill gaps
        Acts::Experimental::detail::BlueprintHelper::fillGaps(*detectorBpr, false);

        auto detectorBuilder =
            std::make_shared<Acts::Experimental::CuboidalContainerBuilder>(
                *detectorBpr, Acts::Logging::VERBOSE);

        auto idGenCfg = LUXEGeometryIdGenerator::Config{
            false, 0u, true, false, gOpt};

        // Detector builder
        Acts::Experimental::DetectorBuilder::Config dCfg;
        dCfg.auxiliary = "LUXE detector builder";
        dCfg.name = "LUXE detector from blueprint";
        dCfg.builder = detectorBuilder;
        dCfg.geoIdGenerator = std::make_shared<LUXEGeometryIdGenerator>(
            idGenCfg,
            Acts::getDefaultLogger("GeoIdGenerator",
                Acts::Logging::VERBOSE));


        auto detector = Acts::Experimental::DetectorBuilder(dCfg).construct(gctx);

        return detector;
}

} // namespace LUXEGeometry
