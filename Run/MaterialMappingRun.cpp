#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/BinnedSurfaceMaterialAccumulater.hpp"
#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Material/MaterialMapper.hpp"

#include <cstddef>

#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryOptions.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/JsonMaterialWriter.hpp"
#include "TrackingPipeline/Io/RootMaterialTrackReader.hpp"
#include "TrackingPipeline/Io/RootMaterialTrackWriter.hpp"
#include "TrackingPipeline/Material/EmptyBinnedMaterialDecorator.hpp"
#include "TrackingPipeline/Material/MaterialMappingAlgorithm.hpp"
#include "toml++/toml.hpp"

using namespace Acts::UnitLiterals;

std::unique_ptr<const E320::GeometryOptions> E320::GeometryOptions::m_instance =
    nullptr;

int main() {
  // Geometry constraints instance
  const auto& goInst = *E320::GeometryOptions::instance();

  // Load configuration
  const std::string pathToCfg =
      "/home/romanurmanov/work/TrackingPipeline/ActsTrackingPipeline/conf/"
      "runs/"
      "MaterialMappingRunConfig.toml";
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
  Acts::BinUtility sensitiveMaterialBinning = goInst.chipMaterialBinningX;
  sensitiveMaterialBinning += goInst.chipMaterialBinningY;

  EmptyBinnedMaterialDecorator::Config materialDecoratorCfg;
  for (const auto& chipPars : goInst.tcParameters) {
    Acts::GeometryIdentifier geoId;
    geoId.setSensitive(chipPars.geoId);
    materialDecoratorCfg.surfaceBinnings.insert(
        {geoId, sensitiveMaterialBinning});
  }

  Acts::BinUtility passiveMaterialBinning = goInst.pdcWindowMaterialBinningX;
  passiveMaterialBinning += goInst.pdcWindowMaterialBinningY;
  Acts::GeometryIdentifier pdcWindowGeoId;
  pdcWindowGeoId.setPassive(goInst.pdcWindowParameters.geoId);
  materialDecoratorCfg.surfaceBinnings.insert(
      {pdcWindowGeoId, passiveMaterialBinning});

  auto materialDecorator =
      std::make_shared<EmptyBinnedMaterialDecorator>(materialDecoratorCfg);

  // Construct detector
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
  // Material mapping setup

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

  // Accumulation setup
  Acts::BinnedSurfaceMaterialAccumulater::Config materialAccumulaterCfg;
  materialAccumulaterCfg.emptyBinCorrection =
      getEntryBool("BinnedSurfaceMaterialAccumulater", "emptyBinCorrection");
  materialAccumulaterCfg.materialSurfaces = materialSurfaces;
  auto materialAccumulater =
      std::make_shared<Acts::BinnedSurfaceMaterialAccumulater>(
          materialAccumulaterCfg,
          Acts::getDefaultLogger("BinnedSurfaceMaterialAccumulater", logLevel));

  // Mapper setup
  Acts::MaterialMapper::Config materialMapperCfg;
  materialMapperCfg.assignmentFinder = materialAssinger;
  materialMapperCfg.surfaceMaterialAccumulater = materialAccumulater;
  auto materialMapper = std::make_shared<Acts::MaterialMapper>(
      materialMapperCfg, Acts::getDefaultLogger("MaterialMapper", logLevel));

  // Initialize material assignment vetos
  double maxAssignmentDownstreamDistance =
      getEntryDouble("MaterialMapper", "maxAssignmentDownstreamDistance");
  double maxAssignmentUpstreamDistance =
      getEntryDouble("MaterialMapper", "maxAssignmentUpstreamDistance");
  Acts::MaterialInteractionAssignment::LocalVeto distanceVeto =
      [&](const Acts::MaterialInteraction& mInt,
          const Acts::IAssignmentFinder::SurfaceAssignment& surfAssignment)
      -> bool {
    double distance = (surfAssignment.position - mInt.position).norm();

    if (mInt.position(goInst.primaryIdx) >
        surfAssignment.position(goInst.primaryIdx)) {
      return (distance > maxAssignmentDownstreamDistance);
    } else {
      return (distance > maxAssignmentUpstreamDistance);
    }
  };

  std::vector<std::pair<Acts::GeometryIdentifier,
                        Acts::MaterialInteractionAssignment::LocalVeto>>
      localVetos;
  for (const auto* surf : materialSurfaces) {
    localVetos.emplace_back(surf->geometryId(), distanceVeto);
  }

  Acts::MaterialMapper::Options materialMapperOpt;
  materialMapperOpt.assignmentOptions.localVetos = Acts::GeometryHierarchyMap<
      Acts::MaterialInteractionAssignment::LocalVeto>(localVetos);

  // Material map writers
  std::vector<std::shared_ptr<IMaterialWriter>> mapWriters;

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

  auto jsonMaterialWriterCfg = JsonMaterialWriter::Config();
  jsonMaterialWriterCfg.converterCfg = jsonMaterialConverterCfg;
  jsonMaterialWriterCfg.filePath =
      getEntryStr("JsonMaterialWriter", "filePath");
  auto jsonMaterialWriter =
      std::make_shared<JsonMaterialWriter>(jsonMaterialWriterCfg, logLevel);

  mapWriters.push_back(jsonMaterialWriter);

  // Mapping Algorithm
  MaterialMappingAlgorithm::Config materialMappingAlgorithmCfg;
  materialMappingAlgorithmCfg.materialMapper = materialMapper;
  materialMappingAlgorithmCfg.materialMapperOptions = materialMapperOpt;
  materialMappingAlgorithmCfg.materiaMaplWriters = mapWriters;
  materialMappingAlgorithmCfg.inputMaterialTracks =
      getEntryStr("MaterialMappingAlgorithm", "inputMaterialTracks");
  materialMappingAlgorithmCfg.mappedMaterialTracks =
      getEntryStr("MaterialMappingAlgorithm", "mappedMaterialTracks");
  materialMappingAlgorithmCfg.unmappedMaterialTracks =
      getEntryStr("MaterialMappingAlgorithm", "unmappedMaterialTracks");
  auto materialMappingAlgorithm = std::make_shared<MaterialMappingAlgorithm>(
      materialMappingAlgorithmCfg, logLevel);

  sequencer.addAlgorithm(materialMappingAlgorithm);

  // Add the mapped material tracks writer
  RootMaterialTrackWriter::Config mappedMaterialTrackWriterCfg;
  mappedMaterialTrackWriterCfg.inputMaterialTracks =
      getEntryStr("RootMaterialTrackWriter", "inputMaterialTracksMapped");
  mappedMaterialTrackWriterCfg.filePath =
      getEntryStr("RootMaterialTrackWriter", "filePathMapped");
  mappedMaterialTrackWriterCfg.treeName =
      getEntryStr("RootMaterialTrackWriter", "treeNameMapped");

  auto mappedMaterialTrackWriter = std::make_shared<RootMaterialTrackWriter>(
      mappedMaterialTrackWriterCfg, logLevel);

  sequencer.addWriter(mappedMaterialTrackWriter);

  // Add the unmapped material tracks writer
  RootMaterialTrackWriter::Config unmappedMaterialTrackWriterCfg;
  unmappedMaterialTrackWriterCfg.inputMaterialTracks =
      getEntryStr("RootMaterialTrackWriter", "inputMaterialTracksUnmapped");
  unmappedMaterialTrackWriterCfg.filePath =
      getEntryStr("RootMaterialTrackWriter", "filePathUnmapped");
  unmappedMaterialTrackWriterCfg.treeName =
      getEntryStr("RootMaterialTrackWriter", "treeNameUnmapped");

  auto unmappedMaterialTrackWriter = std::make_shared<RootMaterialTrackWriter>(
      unmappedMaterialTrackWriterCfg, logLevel);

  sequencer.addWriter(unmappedMaterialTrackWriter);

  return sequencer.run();
}
