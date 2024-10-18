#include "TrackingPipeline/Geometry/LUXEGeometry.hpp"
#include "TrackingPipeline/Geometry/LUXEGeometryIdGenerator.hpp"
#include "TrackingPipeline/Geometry/LayerBuilderConstruction.hpp"
#include "TrackingPipeline/Material/NoMaterialDecorator.hpp"

#include "Acts/Detector/detail/BlueprintHelper.hpp"
#include "Acts/Detector/CuboidalContainerBuilder.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"
#include "Acts/Plugins/Json/JsonMaterialDecorator.hpp"

#include <vector>

namespace LUXEGeometry {

// TODO: Add the vacuum exit window
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
        std::vector<Acts::BinningValue> trackerBins = {Acts::BinningValue::binZ};

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
            auto layerBuilder =
                makeLayerBuilder<1>(
                    world, gOpt.g4ToWorld, names, 
                    {zBounds}, {Acts::BinningValue::binY});

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

        auto zBounds = std::make_tuple(
            gOpt.dipoleTranslation.z() - gOpt.dipoleBounds[2],
            gOpt.dipoleTranslation.z() + gOpt.dipoleBounds[2]);

        auto layerBuilder =
            makeLayerBuilder<1>(
                world, gOpt.g4ToWorld, names, 
                {zBounds}, {Acts::BinningValue::binZ});

        Acts::Transform3 dipoleTransform = Acts::Transform3::Identity();
        dipoleTransform.rotate(
            gOpt.actsToWorld.rotation().inverse());
        dipoleTransform.translate(gOpt.dipoleTranslation);

        auto dipoleNode = 
            std::make_unique<Acts::Experimental::Blueprint::Node>(
            "Dipole", dipoleTransform,
            Acts::VolumeBounds::eCuboid,
            gOpt.dipoleBounds, layerBuilder);

        trackerBP->add(std::move(dipoleNode));

        return trackerBP;
};

std::shared_ptr<const Acts::Experimental::Detector>
    buildLUXEDetector(
        const std::unique_ptr<Acts::Experimental::Blueprint::Node> 
            detectorBpr,
        const Acts::GeometryContext& gctx,
        const LUXEGeometry::GeometryOptions& gOpt,
        const std::vector<Acts::GeometryIdentifier>& materialVetos) {
            // Complete and fill gaps
            Acts::Experimental::detail::BlueprintHelper::fillGaps(*detectorBpr, false);
    
            auto detectorBuilder =
                std::make_shared<Acts::Experimental::CuboidalContainerBuilder>(
                    *detectorBpr, Acts::Logging::VERBOSE);
    
            auto idGenCfg = LUXEGeometryIdGenerator::Config{
                false, 0u, true, false, gOpt};

            // Initialize the material binning
            Acts::BinUtility materialBinning = gOpt.materialBinningX;
            materialBinning += gOpt.materialBinningY;

            auto mpCfg = NoMaterialDecorator::Config();
            mpCfg.surfaceBinning = materialBinning;
            mpCfg.vetos = materialVetos;

            // Detector builder
            Acts::Experimental::DetectorBuilder::Config dCfg;
            dCfg.auxiliary = "LUXE detector builder";
            dCfg.name = "LUXE detector from blueprint";
            dCfg.builder = detectorBuilder;
            dCfg.materialDecorator = std::make_shared<NoMaterialDecorator>(mpCfg);
            dCfg.geoIdGenerator = std::make_shared<LUXEGeometryIdGenerator>(
                idGenCfg,
                Acts::getDefaultLogger("GeoIdGenerator",
                    Acts::Logging::VERBOSE));
    
            auto detector = Acts::Experimental::DetectorBuilder(dCfg).construct(gctx);
    
            return detector;
}

std::shared_ptr<const Acts::Experimental::Detector>
    buildLUXEDetector(
        const std::unique_ptr<Acts::Experimental::Blueprint::Node> 
            detectorBpr,
        const Acts::GeometryContext& gctx,
        const LUXEGeometry::GeometryOptions& gOpt,
        const std::string jsonMaterialPath,
        const std::vector<Acts::GeometryIdentifier>& materialVetos) {
            // Complete and fill gaps
            Acts::Experimental::detail::BlueprintHelper::fillGaps(*detectorBpr, false);
    
            auto detectorBuilder =
                std::make_shared<Acts::Experimental::CuboidalContainerBuilder>(
                    *detectorBpr, Acts::Logging::VERBOSE);
    
            auto idGenCfg = LUXEGeometryIdGenerator::Config{
                false, 0u, true, false, gOpt};

            // Initialize the material binning
            Acts::BinUtility materialBinning = gOpt.materialBinningX;
            materialBinning += gOpt.materialBinningY;

            // Material provider
            auto jMatMapConverterCfg = 
                Acts::MaterialMapJsonConverter::Config();
            jMatMapConverterCfg.context = gctx;

            auto jMatDecorator = 
                std::make_shared<Acts::JsonMaterialDecorator>(
                    jMatMapConverterCfg, 
                    jsonMaterialPath, 
                    Acts::Logging::VERBOSE);

            // Detector builder
            Acts::Experimental::DetectorBuilder::Config dCfg;
            dCfg.auxiliary = "LUXE detector builder";
            dCfg.name = "LUXE detector from blueprint";
            dCfg.builder = detectorBuilder;
            dCfg.materialDecorator = jMatDecorator;
            dCfg.geoIdGenerator = std::make_shared<LUXEGeometryIdGenerator>(
                idGenCfg,
                Acts::getDefaultLogger("GeoIdGenerator",
                    Acts::Logging::VERBOSE));
    
            auto detector = Acts::Experimental::DetectorBuilder(dCfg, 
                Acts::getDefaultLogger("DetectorBuilder", Acts::Logging::VERBOSE)).construct(gctx);
    
            return detector;
}

} // namespace LUXEGeometry
