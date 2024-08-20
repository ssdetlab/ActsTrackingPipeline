#include "ActsLUXEPipeline/E320Geometry.hpp"
#include "ActsLUXEPipeline/Sequencer.hpp"
#include "ActsLUXEPipeline/ROOTMaterialTrackReader.hpp"
#include "ActsLUXEPipeline/CoreMaterialMapping.hpp"
#include "ActsLUXEPipeline/JsonMaterialWriter.hpp"
#include "ActsLUXEPipeline/ROOTMaterialTrackWriter.hpp"
#include "ActsLUXEPipeline/RandomNumbers.hpp"
#include "ActsLUXEPipeline/MaterialValidation.hpp"

#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Material/BinnedSurfaceMaterialAccumulater.hpp"
#include "Acts/Material/MaterialValidater.hpp"
#include "Acts/Material/MaterialMapper.hpp"

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
        "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/E320Pipeline_gdmls/ettgeom_magnet_pdc_tracker.gdml";
    std::vector<std::string> names{"OPPPSensitive"};

    std::vector<Acts::GeometryIdentifier> materialVeto{};

    std::string materialPath = "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_build/material.json";

    auto trackerBP = 
        E320Geometry::makeBlueprintE320(gdmlPath, names, gOpt);
    auto detector =
        E320Geometry::buildE320Detector(std::move(trackerBP), gctx, gOpt, materialPath, materialVeto);

    for (auto& vol : detector->rootVolumes()) {
        std::cout << "Volume: " << vol->name() << " = " << vol->surfaces().size() << std::endl;
        for (auto& surf : vol->surfaces()) {
            std::cout << "Surface: (" << surf->center(gctx).transpose() << ") = (" << surf->normal(gctx,surf->center(gctx),Acts::Vector3(0,1,0)).transpose() << ")" << std::endl;
            std::cout << "Surface material: " << surf->surfaceMaterial() << std::endl;
        }
    }

    // --------------------------------------------------------------
    // Material validation setup

    // Setup the sequencer
    Sequencer::Config seqCfg;
    seqCfg.events = 100000;
    seqCfg.numThreads = 1;
    Sequencer sequencer(seqCfg);

    // White board for material track reading
    auto whiteBoard = WhiteBoard(Acts::getDefaultLogger("WhiteBoard", logLevel));

    auto rnd = std::make_shared<RandomNumbers>(RandomNumbers::Config());

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

    // Validater setup
    auto matValidaterCfg = Acts::MaterialValidater::Config();
    matValidaterCfg.materialAssigner = materialAssinger;
    auto matValidater = 
        std::make_shared<Acts::MaterialValidater>(
            matValidaterCfg, 
            Acts::getDefaultLogger("MaterialValidater", logLevel));

    // Validation Algorithm
    auto uniform = std::make_shared<UniformVertexGenerator>();
    uniform->mins = Acts::Vector3(-8, 0, 16000);
    uniform->maxs = Acts::Vector3(8, 400, 16000);

    auto matValidationCfg = MaterialValidation::Config();
    matValidationCfg.materialValidater = matValidater;
    matValidationCfg.outputMaterialTracks = "recorded-material-tracks";
    matValidationCfg.ntracks = 1;
    matValidationCfg.etaRange = {4.,4.};
    matValidationCfg.startPosition = uniform;
    matValidationCfg.randomNumberSvc = rnd;
    auto materialValidation = 
       std::make_shared<MaterialValidation>(
            matValidationCfg, 
            logLevel);
    sequencer.addAlgorithm(materialValidation);

    // Add the mapped material tracks writer
    auto mappedMaterialTrackWriterCfg = 
        ROOTMaterialTrackWriter::Config{
            "recorded-material-tracks",
            "recorded-material-tracks.root",
            "RECREATE",
            "recorded-material-tracks"
    };    

    auto mappedMaterialTrackWriter = std::make_shared<
        ROOTMaterialTrackWriter>(
            mappedMaterialTrackWriterCfg, 
            logLevel);

    sequencer.addWriter(mappedMaterialTrackWriter);

    // --------------------------------------------------------------
    // Run all configured algorithms and return the appropriate status.

    return sequencer.run();
}
