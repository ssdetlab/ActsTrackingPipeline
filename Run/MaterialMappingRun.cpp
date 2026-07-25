#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/BinnedSurfaceMaterialAccumulater.hpp"
#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Material/MaterialMapper.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Alignment/detail/AlignmentStoreBuilders.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryOptions.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/JsonMaterialWriter.hpp"
#include "TrackingPipeline/Io/RootMaterialTrackReader.hpp"
#include "TrackingPipeline/Io/RootMaterialTrackWriter.hpp"
#include "TrackingPipeline/Material/CoreMaterialMapping.hpp"
#include "TrackingPipeline/Material/NoMaterialDecorator.hpp"

using namespace Acts::UnitLiterals;

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
  Acts::BinUtility sensitiveMaterialBinning = goInst.chipMaterialBinningX;
  sensitiveMaterialBinning += goInst.chipMaterialBinningY;

  NoMaterialDecorator::Config materialDecoratorCfg;
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

  materialDecoratorCfg.vetos = {};

  auto materialDecorator =
      std::make_shared<NoMaterialDecorator>(materialDecoratorCfg);

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

  auto aStore = detail::makeAlignmentStore(gctx, detector.get());
  AlignmentContext alignCtx(aStore);
  Acts::GeometryContext testCtx{alignCtx};
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive() != 0u) {
        std::cout << "-----------------------------------\n";
        std::cout << "SURFACE " << s->geometryId() << "\n";
        std::cout << "CENTER " << s->center(testCtx).transpose() << " -- "
                  << s->center(Acts::GeometryContext()).transpose() << "\n";
        std::cout << "DELTA "
                  << (s->center(testCtx) - s->center(Acts::GeometryContext()))
                             .transpose() *
                         1e3
                  << "\n";

        std::cout << "NORMAL "
                  << s->normal(testCtx, s->center(testCtx),
                               Acts::Vector3::UnitY())
                         .transpose()
                  << " -- "
                  << s->normal(testCtx, s->center(Acts::GeometryContext()),
                               Acts::Vector3::UnitY())
                         .transpose()
                  << "\n";
        std::cout << "ROTATION \n"
                  << s->transform(testCtx).rotation() << " -- \n"
                  << "\n"
                  << s->transform(Acts::GeometryContext()).rotation() << "\n";

        std::cout << "EXTENT "
                  << s->polyhedronRepresentation(testCtx, 1000).extent()
                  << "\n -- \n"
                  << s->polyhedronRepresentation(Acts::GeometryContext(), 1000)
                         .extent()
                  << "\n";
      }
    }
  }
  gctx = Acts::GeometryContext{alignCtx};

  // --------------------------------------------------------------
  // The magnetic field setup

  auto field = E320::buildMagField(gctx);

  // --------------------------------------------------------------
  // Material mapping setup

  // Setup the sequencer
  Sequencer::Config seqCfg;
  // seqCfg.events = 10000;
  seqCfg.numThreads = 1;
  seqCfg.trackFpes = false;
  Sequencer sequencer(seqCfg);

  // White board for material track reading
  auto whiteBoard = WhiteBoard(Acts::getDefaultLogger("WhiteBoard", logLevel));

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
  materialTrackReaderCfg.outputMaterialTracks = "material-tracks",
  materialTrackReaderCfg.treeName = "material-tracks",
  materialTrackReaderCfg.fileList = {
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_material/"
      "Uniform_DirectZ_TrackerOnly_256x128_1e6/"
      "geant4_material_tracks_recording.root"};
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

  // Accumulation setup : Binned surface material accumulater
  Acts::BinnedSurfaceMaterialAccumulater::Config materialAccumulaterCfg;
  materialAccumulaterCfg.emptyBinCorrection = true;
  materialAccumulaterCfg.materialSurfaces = surfaces;
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

  Acts::MaterialMapper::Options materialMapperOpt;
  materialMapperOpt.assignmentOptions.maxDistance = 10_mm;

  // Material map writers
  std::vector<std::shared_ptr<IMaterialWriter>> mapWriters;

  Acts::MaterialMapJsonConverter::Config jsonMaterialConverterCfg;
  jsonMaterialConverterCfg.context = gctx;
  jsonMaterialConverterCfg.processSensitives = true;
  jsonMaterialConverterCfg.processApproaches = true;
  jsonMaterialConverterCfg.processRepresenting = true;
  jsonMaterialConverterCfg.processBoundaries = true;
  jsonMaterialConverterCfg.processVolumes = true;
  jsonMaterialConverterCfg.processDenseVolumes = false;
  jsonMaterialConverterCfg.processNonMaterial = false;

  auto jsonMaterialWriterCfg = JsonMaterialWriter::Config();
  jsonMaterialWriterCfg.converterCfg = jsonMaterialConverterCfg;
  jsonMaterialWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_material/"
      "material.json";
  auto jsonMaterialWriter =
      std::make_shared<JsonMaterialWriter>(jsonMaterialWriterCfg, logLevel);

  mapWriters.push_back(jsonMaterialWriter);

  // Mapping Algorithm
  CoreMaterialMapping::Config coreMaterialMappingCfg;
  coreMaterialMappingCfg.materialMapper = materialMapper;
  coreMaterialMappingCfg.materialMapperOptions = materialMapperOpt;
  coreMaterialMappingCfg.inputMaterialTracks = "material-tracks";
  coreMaterialMappingCfg.mappedMaterialTracks = "mapped-material-tracks";
  coreMaterialMappingCfg.unmappedMaterialTracks = "unmapped-material-tracks";
  coreMaterialMappingCfg.materiaMaplWriters = mapWriters;
  auto coreMaterialMapping =
      std::make_shared<CoreMaterialMapping>(coreMaterialMappingCfg, logLevel);

  sequencer.addAlgorithm(coreMaterialMapping);

  // Add the mapped material tracks writer
  RootMaterialTrackWriter::Config mappedMaterialTrackWriterCfg;
  mappedMaterialTrackWriterCfg.inputMaterialTracks = "mapped-material-tracks";
  mappedMaterialTrackWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_material/"
      "mapped-material-tracks.root";
  mappedMaterialTrackWriterCfg.fileMode = "RECREATE";
  mappedMaterialTrackWriterCfg.treeName = "mapped-material-tracks";
  mappedMaterialTrackWriterCfg.recalculateTotals = false;
  mappedMaterialTrackWriterCfg.prePostStep = false;
  mappedMaterialTrackWriterCfg.storeSurface = true;
  mappedMaterialTrackWriterCfg.storeVolume = false;
  mappedMaterialTrackWriterCfg.collapseInteractions = false;

  auto mappedMaterialTrackWriter = std::make_shared<RootMaterialTrackWriter>(
      mappedMaterialTrackWriterCfg, logLevel);

  sequencer.addWriter(mappedMaterialTrackWriter);

  // Add the unmapped material tracks writer
  RootMaterialTrackWriter::Config unmappedMaterialTrackWriterCfg;
  unmappedMaterialTrackWriterCfg.inputMaterialTracks =
      "unmapped-material-tracks";
  unmappedMaterialTrackWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_material/"
      "unmapped-material-tracks.root";
  unmappedMaterialTrackWriterCfg.fileMode = "RECREATE";
  unmappedMaterialTrackWriterCfg.treeName = "unmapped-material-tracks";
  unmappedMaterialTrackWriterCfg.recalculateTotals = false;
  unmappedMaterialTrackWriterCfg.prePostStep = false;
  unmappedMaterialTrackWriterCfg.storeSurface = true;
  unmappedMaterialTrackWriterCfg.storeVolume = false;
  unmappedMaterialTrackWriterCfg.collapseInteractions = false;

  auto unmappedMaterialTrackWriter = std::make_shared<RootMaterialTrackWriter>(
      unmappedMaterialTrackWriterCfg, logLevel);

  sequencer.addWriter(unmappedMaterialTrackWriter);

  // --------------------------------------------------------------
  // Run all configured algorithms and return the appropriate status.

  return sequencer.run();
}
