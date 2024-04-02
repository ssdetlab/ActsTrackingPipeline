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
        
        std::size_t numLayers = gOpt.layerZPositions.size();

        // Here binning is dont in the unrotated frame
        // Have to fix the consistency inside Acts
        std::vector<Acts::BinningValue> trackerBins = {Acts::binZ};

        Acts::Transform3 trackerTransform = Acts::Transform3::Identity();
        trackerTransform.rotate(
            gOpt.actsToWorld.rotation().inverse());
        trackerTransform.translate(gOpt.trackerTranslation);

        // Create the tracker node of the blueprint
        auto trackerBP = 
            std::make_unique<Acts::Experimental::Blueprint::Node>(
            "Tracker", trackerTransform,
            Acts::VolumeBounds::eCuboid,
            gOpt.trackerBounds, trackerBins);

        // Iterate over the layers and create 
        // the child nodes
        for (std::size_t i = 0; i < numLayers; i++) {
            
            // Layer bounds
            auto zBounds = std::make_tuple(
                gOpt.layerZPositions[i] - gOpt.layerBounds[2],
                gOpt.layerZPositions[i] + gOpt.layerBounds[2]);

            // As the volumes are already rotated, 
            // the selection has to happen along the y-axis
            auto layerBuilder = makeLayerBuilder<1>(
                world, gOpt.g4ToWorld, names, 
                {zBounds}, {Acts::binY});

            // Convention is that the transformations
            // are with respect to the global frame
            Acts::Vector3 layerTranslation =
                Acts::Vector3(
                    gOpt.armTranslation.x(), 
                    gOpt.armTranslation.y(),
                    gOpt.layerZPositions[i]);

            layerTranslation = 
                gOpt.actsToWorld.rotation().inverse() * layerTranslation; 

            Acts::Transform3 layerTransform = Acts::Transform3::Identity();
            layerTransform.rotate(
                gOpt.actsToWorld.rotation().inverse());
            layerTransform.pretranslate(layerTranslation);

            auto layerNode = 
                std::make_unique<Acts::Experimental::Blueprint::Node>(
                    "layer" + std::to_string(i), layerTransform, 
                    Acts::VolumeBounds::eCuboid, gOpt.layerBounds, 
                    layerBuilder);
    
            trackerBP->add(std::move(layerNode));
        }

        Acts::Transform3 dipoleTransform = Acts::Transform3::Identity();
        dipoleTransform.rotate(
            gOpt.actsToWorld.rotation().inverse());
        dipoleTransform.translate(gOpt.dipoleTranslation);

        auto dipoleBP = 
            std::make_unique<Acts::Experimental::Blueprint::Node>(
            "Dipole", dipoleTransform,
            Acts::VolumeBounds::eCuboid,
            gOpt.dipoleBounds, trackerBins);

        trackerBP->add(std::move(dipoleBP));

        return trackerBP;
};

std::shared_ptr<const Acts::Experimental::Detector>
    buildLUXEDetector(
        const std::unique_ptr<Acts::Experimental::Blueprint::Node> 
            detectorBpr,
        const Acts::GeometryContext& gctx,
        const LUXEGeometry::GeometryOptions& gOpt) {
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
