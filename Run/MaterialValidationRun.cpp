#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Material/MaterialValidater.hpp"
#include "Acts/Plugins/Json/JsonMaterialDecorator.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryOptions.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/RootMaterialTrackReader.hpp"
#include "TrackingPipeline/Io/RootMaterialTrackWriter.hpp"
#include "TrackingPipeline/Material/MaterialValidationAlgorithm.hpp"
#include "toml++/toml.hpp"

std::unique_ptr<const E320::GeometryOptions> E320::GeometryOptions::m_instance =
    nullptr;

int main() {
  // Geometry constraints instance
  const auto& goInst = *E320::GeometryOptions::instance();

  // Load configuration
  const std::string pathToCfg =
      "/home/romanurmanov/work/TrackingPipeline/ActsTrackingPipeline/conf/"
      "runs/"
      "MaterialValidationRun.toml";
  auto runCfg = toml::parse_file(pathToCfg);

  auto getEntryDouble = [&runCfg](const std::string& section,
                                  const std::string& subsection) {
    return runCfg[section][subsection].value<double>().value();
  };
  auto getEntrySizeT = [&runCfg](const std::string& section,
                                 const std::string& subsection) {
    return runCfg[section][subsection].value<std::size_t>().value();
  };
  auto getEntryInt = [&runCfg](const std::string& section,
                               const std::string& subsection) {
    return runCfg[section][subsection].value<int>().value();
  };
  auto getEntryBool = [&runCfg](const std::string& section,
                                const std::string& subsection) {
    return runCfg[section][subsection].value<bool>().value();
  };
  auto getEntryStr = [&runCfg](const std::string& section,
                               const std::string& subsection) {
    return runCfg[section][subsection].value<std::string>().value();
  };

  // Set the log level
  Acts::Logging::Level logLevel =
      Acts::Logging::Level(getEntrySizeT("General", "logLevel"));

  // Initialize contexts
  Acts::GeometryContext gctx;

  // --------------------------------------------------------------
  // Detector setup

  // Material decorator
  Acts::MaterialMapJsonConverter::Config jsonMaterialConverterCfg;
  jsonMaterialConverterCfg.context = gctx;
  jsonMaterialConverterCfg.processSensitives =
      getEntryBool("MaterialMapJsonConverter", "processSensitives");
  jsonMaterialConverterCfg.processApproaches =
      getEntryBool("MaterialMapJsonConverter", "processApproaches");
  jsonMaterialConverterCfg.processRepresenting =
      getEntryBool("MaterialMapJsonConverter", "processRepresenting");
  jsonMaterialConverterCfg.processBoundaries =
      getEntryBool("MaterialMapJsonConverter", "processBoundaries");
  jsonMaterialConverterCfg.processVolumes =
      getEntryBool("MaterialMapJsonConverter", "processVolumes");
  jsonMaterialConverterCfg.processDenseVolumes =
      getEntryBool("MaterialMapJsonConverter", "processDenseVolumes");
  jsonMaterialConverterCfg.processNonMaterial =
      getEntryBool("MaterialMapJsonConverter", "processNonMaterial");

  auto materialDecorator = std::make_shared<Acts::JsonMaterialDecorator>(
      jsonMaterialConverterCfg,
      getEntryStr("JsonMaterialDecorator", "materialPath"), logLevel);

  // Construct detector
  auto detector = getEntryBool("Geometry", "materialDecorator")
                      ? E320::buildDetector(gctx, materialDecorator)
                      : E320::buildDetector(gctx, nullptr);

  for (const auto& vol : detector->volumes()) {
    std::cout << "------------------------------------------\n";
    std::cout << vol->name() << "\n";
    std::cout << vol->extent(gctx);
    std::cout << "Surfaces:\n";
    for (const auto& surf : vol->surfaces()) {
      std::cout << surf->geometryId() << "\n";
      std::cout << surf->center(gctx) << "\n";
      std::cout << surf->polyhedronRepresentation(gctx, 1000).extent() << "\n";
    }
  }

  // --------------------------------------------------------------
  // The magnetic field setup

  auto field = E320::buildMagField(gctx);

  // --------------------------------------------------------------
  // Material validation setup

  // Setup the sequencer
  Sequencer::Config seqCfg;
  seqCfg.skip = getEntrySizeT("Sequencer", "skip");
  seqCfg.events = getEntrySizeT("Sequencer", "events");
  seqCfg.numThreads = 1;
  seqCfg.trackFpes = false;
  seqCfg.logLevel = logLevel;
  Sequencer sequencer(seqCfg);

  // Add the material track reader
  Acts::Transform3 toWorldTransform = Acts::Transform3::Identity();

  Acts::RotationMatrix3 surfToWorldRotationLong =
      Acts::AngleAxis3(M_PI_2, goInst.longDir).toRotationMatrix();
  toWorldTransform.rotate(surfToWorldRotationLong);

  RootMaterialTrackReader::Config materialTrackReaderCfg;
  materialTrackReaderCfg.toWorldTransform = toWorldTransform;
  materialTrackReaderCfg.outputMaterialTracks =
      getEntryStr("RootMaterialTrackReader", "outputMaterialTracks");
  materialTrackReaderCfg.treeName =
      getEntryStr("RootMaterialTrackReader", "treeName");
  materialTrackReaderCfg.readCachedSurfaceInformation =
      getEntryBool("RootMaterialTrackReader", "readCachedSurfaceInformation");
  materialTrackReaderCfg.filePaths = {
      getEntryStr("RootMaterialTrackReader", "filePath")};

  auto materialTrackReader = std::make_shared<RootMaterialTrackReader>(
      materialTrackReaderCfg, logLevel);

  sequencer.addReader(materialTrackReader);

  // Assignment setup
  Acts::IntersectionMaterialAssigner::Config materialAssingerCfg;
  std::vector<const Acts::Surface*> materialSurfaces;
  for (auto& vol : detector->rootVolumes()) {
    for (auto& surf : vol->surfaces()) {
      const auto& geoId = surf->geometryId();
      if (geoId.sensitive() >= goInst.tcParameters.front().geoId &&
              geoId.sensitive() <= goInst.tcParameters.back().geoId ||
          geoId.passive() == goInst.pdcWindowParameters.geoId) {
        materialSurfaces.push_back(surf);
        std::cout << "MATERIAL Surface: (" << surf->center(gctx).transpose()
                  << ") = ("
                  << surf->normal(gctx, surf->center(gctx),
                                  Acts::Vector3(0, 1, 0))
                         .transpose()
                  << ")" << std::endl;
      }
    }
  }

  materialAssingerCfg.surfaces = materialSurfaces;
  auto materialAssinger = std::make_shared<Acts::IntersectionMaterialAssigner>(
      materialAssingerCfg,
      Acts::getDefaultLogger("IntersectionMaterialAssigner", logLevel));

  // Validater setup
  Acts::MaterialValidater::Config matValidaterCfg;
  matValidaterCfg.materialAssigner = materialAssinger;

  auto matValidater = std::make_shared<Acts::MaterialValidater>(
      matValidaterCfg, Acts::getDefaultLogger("MaterialValidater", logLevel));

  // Validation setup
  MaterialValidationAlgorithm::Config materialValidationAlgorithmCfg;
  materialValidationAlgorithmCfg.materialValidater = matValidater;
  materialValidationAlgorithmCfg.inputMaterialTracks =
      getEntryStr("MaterialValidationAlgorithm", "inputMaterialTracks");
  materialValidationAlgorithmCfg.outputMaterialTracks =
      getEntryStr("MaterialValidationAlgorithm", "outputMaterialTracks");

  auto materialValidationAlgorithm =
      std::make_shared<MaterialValidationAlgorithm>(
          materialValidationAlgorithmCfg, logLevel);
  sequencer.addAlgorithm(materialValidationAlgorithm);

  // Add the recorded material tracks writer
  RootMaterialTrackWriter::Config recordedMaterialTrackWriterCfg;
  recordedMaterialTrackWriterCfg.inputMaterialTracks =
      getEntryStr("RootMaterialTrackWriter", "inputMaterialTracks");
  recordedMaterialTrackWriterCfg.filePath =
      getEntryStr("RootMaterialTrackWriter", "filePath");
  recordedMaterialTrackWriterCfg.treeName =
      getEntryStr("RootMaterialTrackWriter", "treeName");

  auto recordedMaterialTrackWriter = std::make_shared<RootMaterialTrackWriter>(
      recordedMaterialTrackWriterCfg, logLevel);

  sequencer.addWriter(recordedMaterialTrackWriter);

  // --------------------------------------------------------------
  // Run all configured algorithms and return the appropriate status.

  return sequencer.run();
}
