#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/Sequencer.hpp"
#include "ActsLUXEPipeline/ROOTMaterialTrackReader.hpp"
#include "ActsLUXEPipeline/CoreMaterialMapping.hpp"
#include "ActsLUXEPipeline/JsonMaterialWriter.hpp"
#include "ActsLUXEPipeline/ROOTMaterialTrackWriter.hpp"

#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Material/BinnedSurfaceMaterialAccumulater.hpp"
#include "Acts/Material/MaterialMapper.hpp"

int main() {
    // Set the log level
    Acts::Logging::Level logLevel = Acts::Logging::INFO;

    // Dummy context and options
    Acts::GeometryContext gctx;
    Acts::MagneticFieldContext mctx;
    Acts::CalibrationContext cctx;
    LUXEGeometry::GeometryOptions gOpt;

    // --------------------------------------------------------------
    // LUXE detector setup

    // Set the path to the gdml file
    // and the names of the volumes to be converted
    std::string gdmlPath = 
        "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_gdmls/lxgeomdump_ip_tracker_positron.gdml";
    std::vector<std::string> names{"OPPPSensitive", "VCWindowPanel"};

    std::vector<Acts::GeometryIdentifier> materialVeto{};

    auto trackerBP = 
        LUXEGeometry::makeBlueprintLUXE(gdmlPath, names, gOpt);
    auto detector =
        LUXEGeometry::buildLUXEDetector(std::move(trackerBP), gctx, gOpt, materialVeto);

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
    // seqCfg.events = 10;
    seqCfg.numThreads = 1;
    Sequencer sequencer(seqCfg);

    // White board for material track reading
    auto whiteBoard = WhiteBoard(Acts::getDefaultLogger("WhiteBoard", logLevel));

    // Add the material track reader
    auto materialTrackReaderCfg = ROOTMaterialTrackReader::Config{
        "material-tracks",
        "material-tracks",
        {"/home/romanurmanov/lab/LUXE/reports/acts_telescope_geo_4/material_files/uniform/geant4_material_tracks.root"}
    };

    auto materialTrackReader = std::make_shared<ROOTMaterialTrackReader>(
        materialTrackReaderCfg,
        logLevel);

    sequencer.addReader(materialTrackReader);

    // Assignment setup : Intersection assigner
    auto materialAssingerCfg = Acts::IntersectionMaterialAssigner::Config();
    std::vector<const Acts::Surface*> surfaces;
    for (auto& vol : detector->rootVolumes()) {
        for (auto& surf : vol->surfaces()) {
            surfaces.push_back(surf);
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
        ROOTMaterialTrackWriter::Config{
            "mapped-material-tracks",
            "mapped-material-tracks.root",
            "RECREATE",
            "mapped-material-tracks"
    };    

    auto mappedMaterialTrackWriter = std::make_shared<
        ROOTMaterialTrackWriter>(
            mappedMaterialTrackWriterCfg, 
            logLevel);

    sequencer.addWriter(mappedMaterialTrackWriter);

    // Add the unmapped material tracks writer
    auto unmappedMaterialTrackWriterCfg = 
        ROOTMaterialTrackWriter::Config{
            "unmapped-material-tracks",
            "unmapped-material-tracks.root",
            "RECREATE",
            "unmapped-material-tracks"
    };

    auto unmappedMaterialTrackWriter = std::make_shared<
        ROOTMaterialTrackWriter>(
            unmappedMaterialTrackWriterCfg, 
            logLevel);

    sequencer.addWriter(unmappedMaterialTrackWriter);

    // --------------------------------------------------------------
    // Run all configured algorithms and return the appropriate status.

    return sequencer.run();
}
