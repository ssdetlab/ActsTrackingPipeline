#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Material/MaterialValidater.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/RootMaterialTrackReader.hpp"
#include "TrackingPipeline/Io/RootMaterialTrackWriter.hpp"
#include "TrackingPipeline/Material/MaterialValidation.hpp"

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
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_gdmls/"
      "ettgeom_magnet_pdc_tracker.gdml";
  std::vector<std::string> names{"OPPPSensitive", "DetChamberWindow"};

  // Veto PDC window material mapping
  // to preserve homogeneous material
  // from Geant4
  Acts::GeometryIdentifier pdcWindowId;
  pdcWindowId.setApproach(1);
  std::vector<Acts::GeometryIdentifier> materialVeto{pdcWindowId};

  std::string materialPath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_material/"
      "Uniform_DirectZ_TrackerOnly_256x128_1M/material.json";

  auto trackerBP = E320Geometry::makeBlueprintE320(gdmlPath, names, gOpt);
  auto detector = E320Geometry::buildE320Detector(
      std::move(trackerBP), gctx, gOpt, materialPath, materialVeto);

  for (auto& vol : detector->rootVolumes()) {
    std::cout << "Volume: " << vol->name() << " = " << vol->surfaces().size()
              << std::endl;
    for (auto& surf : vol->surfaces()) {
      std::cout << "Surface: (" << surf->center(gctx).transpose() << ") = ("
                << surf->normal(gctx, surf->center(gctx),
                                Acts::Vector3(0, 1, 0))
                       .transpose()
                << ")" << surf->geometryId() << std::endl;
      std::cout << "Surface material: " << surf->surfaceMaterial() << std::endl;
    }
  }

  // --------------------------------------------------------------
  // Material validation setup

  // Setup the sequencer
  Sequencer::Config seqCfg;
  // seqCfg.events = 1;
  seqCfg.numThreads = 1;
  Sequencer sequencer(seqCfg);

  // Add the material track reader
  auto materialTrackReaderCfg = RootMaterialTrackReader::Config{
      "materialTracks",
      "material-tracks",
      {"/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_material/"
       "Uniform_DirectZ_TrackerOnly_256x128_1M/"
       "geant4_material_tracks_validation.root"}};

  auto materialTrackReader = std::make_shared<RootMaterialTrackReader>(
      materialTrackReaderCfg, logLevel);

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
  auto materialAssinger = std::make_shared<Acts::IntersectionMaterialAssigner>(
      materialAssingerCfg,
      Acts::getDefaultLogger("IntersectionMaterialAssigner", logLevel));

  // Validater setup
  auto matValidaterCfg = Acts::MaterialValidater::Config();
  matValidaterCfg.materialAssigner = materialAssinger;
  auto matValidater = std::make_shared<Acts::MaterialValidater>(
      matValidaterCfg, Acts::getDefaultLogger("MaterialValidater", logLevel));

  // Validation Algorithm
  auto matValidationCfg = MaterialValidation::Config();
  matValidationCfg.materialValidater = matValidater;
  matValidationCfg.inputMaterialTracks = "materialTracks";
  matValidationCfg.outputMaterialTracks = "recordedMaterialTracks";
  auto materialValidation =
      std::make_shared<MaterialValidation>(matValidationCfg, logLevel);
  sequencer.addAlgorithm(materialValidation);

  // Add the mapped material tracks writer
  auto mappedMaterialTrackWriterCfg = RootMaterialTrackWriter::Config{
      "recordedMaterialTracks", "recorded-material-tracks.root", "RECREATE",
      "recorded-material-tracks"};

  auto mappedMaterialTrackWriter = std::make_shared<RootMaterialTrackWriter>(
      mappedMaterialTrackWriterCfg, logLevel);

  sequencer.addWriter(mappedMaterialTrackWriter);

  // --------------------------------------------------------------
  // Run all configured algorithms and return the appropriate status.

  return sequencer.run();
}
