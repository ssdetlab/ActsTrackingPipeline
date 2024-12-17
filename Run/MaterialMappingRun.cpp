#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/RootMaterialTrackReader.hpp"
#include "TrackingPipeline/Io/JsonMaterialWriter.hpp"
#include "TrackingPipeline/Io/RootMaterialTrackWriter.hpp"
#include "TrackingPipeline/Material/CoreMaterialMapping.hpp"

#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Material/BinnedSurfaceMaterialAccumulater.hpp"
#include "Acts/Material/MaterialMapper.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

int main() {
    // Set the log level
    Acts::Logging::Level logLevel = Acts::Logging::INFO;

    // Dummy context and options
    Acts::GeometryContext gctx;
    Acts::MagneticFieldContext mctx;
    Acts::CalibrationContext cctx;
    E320Geometry::GeometryOptions gOpt;

    // --------------------------------------------------------------
    // Detector setup

    // Set the path to the gdml file
    // and the names of the volumes to be converted
    std::string gdmlPath = 
        "/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_gdmls/ettgeom_magnet_pdc_tracker.gdml";
    std::vector<std::string> names{"OPPPSensitive", "DetChamberWindow"};

    std::vector<Acts::GeometryIdentifier> materialVeto{};

    auto trackerBP = 
        E320Geometry::makeBlueprintE320(gdmlPath, names, gOpt);
    auto detector =
        E320Geometry::buildE320Detector(std::move(trackerBP), gctx, gOpt, materialVeto);

    for (auto& vol : detector->rootVolumes()) {
        std::cout << "Volume: " << vol->name() << " = " << vol->surfaces().size() << std::endl;
        for (auto& surf : vol->surfaces()) {
            std::cout << "Surface: (" << surf->center(gctx).transpose() << ") = (" << surf->normal(gctx,surf->center(gctx),Acts::Vector3(0,1,0)).transpose() << ")" << std::endl;
        }
    }

    // --------------------------------------------------------------
    // Material mapping setup

    // Setup the sequencer
    Sequencer::Config seqCfg;
    seqCfg.numThreads = 1;
    seqCfg.trackFpes = false;
    Sequencer sequencer(seqCfg);

    // White board for material track reading
    auto whiteBoard = WhiteBoard(Acts::getDefaultLogger("WhiteBoard", logLevel));

    // Add the material track reader
    auto materialTrackReaderCfg = RootMaterialTrackReader::Config{
        "material-tracks",
        "material-tracks",
        {"/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_material/Uniform_DirectZ_TrackerOnly_256x128_1M/geant4_material_tracks_mapping.root"}
    };

    auto materialTrackReader = std::make_shared<RootMaterialTrackReader>(
        materialTrackReaderCfg,
        logLevel);

    sequencer.addReader(materialTrackReader);

    // Assignment setup : Intersection assigner
    auto materialAssingerCfg = Acts::IntersectionMaterialAssigner::Config();
    std::vector<const Acts::Surface*> surfaces;
    for (auto& vol : detector->rootVolumes()) {
        for (auto& surf : vol->surfaces()) {
            surfaces.push_back(surf);
            std::cout << "MATERIAL Surface: (" << surf->center(gctx).transpose() << ") = (" << surf->normal(gctx,surf->center(gctx),Acts::Vector3(0,1,0)).transpose() << ")" << std::endl;
        }
    }

    materialAssingerCfg.surfaces = surfaces;
    auto materialAssinger = std::make_shared<
        Acts::IntersectionMaterialAssigner>(
            materialAssingerCfg, 
            Acts::getDefaultLogger("IntersectionMaterialAssigner", logLevel));

    // Accumulation setup : Binned surface material accumulater
    auto materialAccumulaterCfg = Acts::BinnedSurfaceMaterialAccumulater::Config();
    materialAccumulaterCfg.materialSurfaces = surfaces;
    auto materialAccumulater = std::make_shared<
        Acts::BinnedSurfaceMaterialAccumulater>(
            materialAccumulaterCfg, 
            Acts::getDefaultLogger("BinnedSurfaceMaterialAccumulater", logLevel));

    // Mapper setup
    auto materialMapperCfg = Acts::MaterialMapper::Config();
    materialMapperCfg.assignmentFinder = materialAssinger;
    materialMapperCfg.surfaceMaterialAccumulater = materialAccumulater;
    auto materialMapper = std::make_shared<
        Acts::MaterialMapper>(
            materialMapperCfg, 
            Acts::getDefaultLogger("MaterialMapper", logLevel));

    // Material map writers
    std::vector<std::shared_ptr<IMaterialWriter>> mapWriters;

    auto jsonMaterialWriterCfg = JsonMaterialWriter::Config();
    auto jsonMaterialWriter = std::make_shared<JsonMaterialWriter>(
        jsonMaterialWriterCfg, 
        logLevel);

    mapWriters.push_back(jsonMaterialWriter);

    // Mapping Algorithm
    auto coreMaterialMappingCfg = CoreMaterialMapping::Config();
    coreMaterialMappingCfg.materialMapper = materialMapper;
    coreMaterialMappingCfg.inputMaterialTracks = "material-tracks";
    coreMaterialMappingCfg.mappedMaterialTracks = "mapped-material-tracks";
    coreMaterialMappingCfg.unmappedMaterialTracks = "unmapped-material-tracks";
    coreMaterialMappingCfg.materiaMaplWriters = mapWriters;
    auto coreMaterialMapping = std::make_shared<
        CoreMaterialMapping>(
            coreMaterialMappingCfg, 
            logLevel);
    
    sequencer.addAlgorithm(coreMaterialMapping);

    // Add the mapped material tracks writer
    auto mappedMaterialTrackWriterCfg = 
        RootMaterialTrackWriter::Config{
            "mapped-material-tracks",
            "mapped-material-tracks.root",
            "RECREATE",
            "mapped-material-tracks"
    };    

    auto mappedMaterialTrackWriter = std::make_shared<
        RootMaterialTrackWriter>(
            mappedMaterialTrackWriterCfg, 
            logLevel);

    sequencer.addWriter(mappedMaterialTrackWriter);

    // Add the unmapped material tracks writer
    auto unmappedMaterialTrackWriterCfg = 
        RootMaterialTrackWriter::Config{
            "unmapped-material-tracks",
            "unmapped-material-tracks.root",
            "RECREATE",
            "unmapped-material-tracks"
    };

    auto unmappedMaterialTrackWriter = std::make_shared<
        RootMaterialTrackWriter>(
            unmappedMaterialTrackWriterCfg, 
            logLevel);

    sequencer.addWriter(unmappedMaterialTrackWriter);

    // --------------------------------------------------------------
    // Run all configured algorithms and return the appropriate status.

    return sequencer.run();
}
