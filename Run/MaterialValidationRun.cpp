#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Material/MaterialValidater.hpp"
#include "Acts/Plugins/Json/JsonMaterialDecorator.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/RootMaterialTrackReader.hpp"
#include "TrackingPipeline/Io/RootMaterialTrackWriter.hpp"
#include "TrackingPipeline/Material/MaterialValidation.hpp"

std::unique_ptr<const E320::GeometryOptions> E320::GeometryOptions::m_instance =
    nullptr;

int main() {
  const auto& goInst = *E320::GeometryOptions::instance();

  // Set the log level
  Acts::Logging::Level logLevel = Acts::Logging::INFO;

  // Dummy context and options
  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext cctx;

  // --------------------------------------------------------------
  // Detector setup

  // Material decorator
  Acts::MaterialMapJsonConverter::Config jsonMaterialConverterCfg;
  jsonMaterialConverterCfg.context = gctx;
  jsonMaterialConverterCfg.processSensitives = true;
  jsonMaterialConverterCfg.processApproaches = true;
  jsonMaterialConverterCfg.processRepresenting = true;
  jsonMaterialConverterCfg.processBoundaries = true;
  jsonMaterialConverterCfg.processVolumes = true;
  jsonMaterialConverterCfg.processDenseVolumes = false;
  jsonMaterialConverterCfg.processNonMaterial = false;

  auto materialDecorator = std::make_shared<Acts::JsonMaterialDecorator>(
      jsonMaterialConverterCfg,
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_material/"
      "Uniform_DirectZ_Tracker_PDCWindow_256x128_1e6/material.json",
      logLevel);

  auto detector = E320::buildDetector(gctx, materialDecorator);

  std::map<Acts::GeometryIdentifier, const Acts::Surface*> surfaceMap;
  for (const auto& vol : detector->volumes()) {
    std::cout << "------------------------------------------\n";
    std::cout << vol->name() << "\n";
    std::cout << vol->extent(gctx);
    std::cout << "Surfaces:\n";
    for (const auto& surf : vol->surfaces()) {
      std::cout << surf->geometryId() << "\n";
      std::cout << surf->center(gctx) << "\n";
      std::cout << surf->polyhedronRepresentation(gctx, 1000).extent() << "\n";
      if (surf->geometryId().sensitive() != 0u) {
        surfaceMap[surf->geometryId()] = surf;
      }
    }
  }

  // --------------------------------------------------------------
  // The magnetic field setup

  auto field = E320::buildMagField(gctx);

  // --------------------------------------------------------------
  // Material validation setup

  // Setup the sequencer
  Sequencer::Config seqCfg;
  // seqCfg.events = 1;
  seqCfg.numThreads = 1;
  Sequencer sequencer(seqCfg);

  // Add the material track reader
  Acts::Transform3 toWorldTransform = Acts::Transform3::Identity();

  Acts::RotationMatrix3 toWorldRotationX =
      Acts::AngleAxis3(goInst.toWorldAngleX, Acts::Vector3::UnitX())
          .toRotationMatrix();
  Acts::RotationMatrix3 toWorldRotationY =
      Acts::AngleAxis3(M_PI_2, Acts::Vector3::UnitY()).toRotationMatrix();
  Acts::RotationMatrix3 toWorldRotationZ =
      Acts::AngleAxis3(0, Acts::Vector3::UnitZ()).toRotationMatrix();

  toWorldTransform.translate(Acts::Vector3::Zero());

  toWorldTransform.rotate(toWorldRotationX);
  toWorldTransform.rotate(toWorldRotationY);
  toWorldTransform.rotate(toWorldRotationZ);

  RootMaterialTrackReader::Config materialTrackReaderCfg;
  materialTrackReaderCfg.outputMaterialTracks = "material-tracks";
  materialTrackReaderCfg.treeName = "material-tracks";
  materialTrackReaderCfg.fileList = {
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_material/"
      "Uniform_DirectZ_Tracker_PDCWindow_256x128_1e6/"
      "geant4-material-tracks-validation.root"};
  materialTrackReaderCfg.toWorldTransform = toWorldTransform;
  materialTrackReaderCfg.readCachedSurfaceInformation = false;

  auto materialTrackReader = std::make_shared<RootMaterialTrackReader>(
      materialTrackReaderCfg, logLevel);

  sequencer.addReader(materialTrackReader);

  // Assignment setup : Intersection assigner
  Acts::IntersectionMaterialAssigner::Config materialAssingerCfg;
  std::vector<const Acts::Surface*> surfaces;
  for (auto& vol : detector->rootVolumes()) {
    for (auto& surf : vol->surfaces()) {
      if (surf->geometryId().sensitive() >= 40) {
        continue;
      }
      surfaces.push_back(surf);
      std::cout << "MATERIAL Surface: (" << surf->center(gctx).transpose()
                << ") = ("
                << surf->normal(gctx, surf->center(gctx),
                                Acts::Vector3(0, 1, 0))
                       .transpose()
                << ")" << std::endl;
    }
  }

  materialAssingerCfg.surfaces = surfaces;
  auto materialAssinger = std::make_shared<Acts::IntersectionMaterialAssigner>(
      materialAssingerCfg,
      Acts::getDefaultLogger("IntersectionMaterialAssigner", logLevel));

  // Validater setup
  Acts::MaterialValidater::Config matValidaterCfg;
  matValidaterCfg.materialAssigner = materialAssinger;

  auto matValidater = std::make_shared<Acts::MaterialValidater>(
      matValidaterCfg, Acts::getDefaultLogger("MaterialValidater", logLevel));

  // Validation Algorithm
  MaterialValidation::Config matValidationCfg;
  matValidationCfg.materialValidater = matValidater;
  matValidationCfg.inputMaterialTracks = "material-tracks";
  matValidationCfg.outputMaterialTracks = "recorded-material-tracks";
  auto materialValidation =
      std::make_shared<MaterialValidation>(matValidationCfg, logLevel);
  sequencer.addAlgorithm(materialValidation);

  // Add the recorded material tracks writer
  RootMaterialTrackWriter::Config recordedMaterialTrackWriterCfg;
  recordedMaterialTrackWriterCfg.inputMaterialTracks =
      "recorded-material-tracks";
  recordedMaterialTrackWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_material/"
      "recorded-material-tracks.root";
  recordedMaterialTrackWriterCfg.fileMode = "RECREATE";
  recordedMaterialTrackWriterCfg.treeName = "recorded-material-tracks";
  recordedMaterialTrackWriterCfg.recalculateTotals = false;
  recordedMaterialTrackWriterCfg.prePostStep = false;
  recordedMaterialTrackWriterCfg.storeSurface = true;
  recordedMaterialTrackWriterCfg.storeVolume = false;
  recordedMaterialTrackWriterCfg.collapseInteractions = false;

  auto recordedMaterialTrackWriter = std::make_shared<RootMaterialTrackWriter>(
      recordedMaterialTrackWriterCfg, logLevel);

  sequencer.addWriter(recordedMaterialTrackWriter);

  // --------------------------------------------------------------
  // Run all configured algorithms and return the appropriate status.

  return sequencer.run();
}
