#include "TrackingPipeline/Geometry/E320Geometry.hpp"

#include "Acts/Detector/CuboidalContainerBuilder.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"
#include "Acts/Detector/detail/BlueprintHelper.hpp"
#include "Acts/Plugins/Json/JsonMaterialDecorator.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <memory>
#include <vector>

#include "TrackingPipeline/Alignment/AlignableDetectorElement.hpp"
#include "TrackingPipeline/Alignment/AlignableDetectorElementBuilder.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/E320GeometryIdGenerator.hpp"
#include "TrackingPipeline/Geometry/LayerBuilderConstruction.hpp"
#include "TrackingPipeline/Material/NoMaterialDecorator.hpp"

namespace E320Geometry {

std::unique_ptr<Acts::Experimental::Blueprint::Node> makeBlueprintE320(
    const std::string& gdmlPath, const std::vector<std::string>& names,
    const E320Geometry::GeometryOptions& gOpt) {
  // Read the gdml file and get the world volume
  G4GDMLParser parser;
  parser.Read(gdmlPath, false);
  auto world = parser.GetWorldVolume();

  std::size_t numLayers = gOpt.layerZPositions.size();

  // Here binning is dont in the unrotated frame
  // Have to fix the consistency inside Acts
  std::vector<Acts::BinningValue> trackerBins = {Acts::BinningValue::binZ};

  Acts::Transform3 trackerTransform = Acts::Transform3::Identity();
  trackerTransform.rotate(gOpt.actsToWorld.rotation().inverse());
  trackerTransform.translate(gOpt.trackerTranslation);

  // Create the tracker node of the blueprint
  auto trackerBP = std::make_unique<Acts::Experimental::Blueprint::Node>(
      "Tracker", trackerTransform, Acts::VolumeBounds::eCuboid,
      gOpt.trackerBounds, trackerBins);

  // Iterate over the layers and create
  // the child nodes
  for (std::size_t i = 0; i < numLayers; i++) {
    // Layer bounds
    auto zBounds =
        std::make_tuple(gOpt.layerZPositions.at(i) - gOpt.layerBounds.at(2),
                        gOpt.layerZPositions.at(i) + gOpt.layerBounds.at(2));

    // As the volumes are already rotated,
    // the selection has to happen along the y-axis
    auto layerBuilder = makeLayerBuilder(world, gOpt.g4ToWorld, names,
                                         {zBounds}, {Acts::BinningValue::binY});

    auto detElementBuilder =
        std::make_shared<AlignableDetectorElementBuilder>();

    // Convention is that the transformations
    // are with respect to the global frame
    Acts::Vector3 layerTranslation =
        Acts::Vector3(gOpt.armTranslation.x(), gOpt.armTranslation.y(),
                      gOpt.layerZPositions.at(i));

    layerTranslation = gOpt.actsToWorld.rotation().inverse() * layerTranslation;

    Acts::Transform3 layerTransform = Acts::Transform3::Identity();
    layerTransform.rotate(gOpt.actsToWorld.rotation().inverse());
    layerTransform.pretranslate(layerTranslation);

    auto layerNode = std::make_unique<Acts::Experimental::Blueprint::Node>(
        "layer" + std::to_string(i), layerTransform,
        Acts::VolumeBounds::eCuboid, gOpt.layerBounds, layerBuilder,
        detElementBuilder);

    trackerBP->add(std::move(layerNode));
  }

  // PDC window
  auto pdcWindowZBounds = std::make_tuple(
      gOpt.pdcWindowTranslation.z() - gOpt.pdcWindowBounds.at(2),
      gOpt.pdcWindowTranslation.z() + gOpt.pdcWindowBounds.at(2));

  auto pdcWindowLayerBuilder =
      makeLayerBuilder(world, gOpt.g4ToWorld, names, {pdcWindowZBounds},
                       {Acts::BinningValue::binY}, true);

  Acts::Transform3 pdcWindowTransform = Acts::Transform3::Identity();
  pdcWindowTransform.rotate(gOpt.actsToWorld.rotation().inverse());
  pdcWindowTransform.translate(gOpt.pdcWindowTranslation);

  auto pdcWindowNode = std::make_unique<Acts::Experimental::Blueprint::Node>(
      "pdcWindow", pdcWindowTransform, Acts::VolumeBounds::eCuboid,
      gOpt.pdcWindowBounds, pdcWindowLayerBuilder);

  trackerBP->add(std::move(pdcWindowNode));

  // Dipole
  auto dipoleZBounds =
      std::make_tuple(gOpt.dipoleTranslation.z() - gOpt.dipoleBounds.at(2),
                      gOpt.dipoleTranslation.z() + gOpt.dipoleBounds.at(2));

  auto dipoleLayerBuilder =
      makeLayerBuilder(world, gOpt.g4ToWorld, names, {dipoleZBounds},
                       {Acts::BinningValue::binZ});

  Acts::Transform3 dipoleTransform = Acts::Transform3::Identity();
  dipoleTransform.rotate(gOpt.actsToWorld.rotation().inverse());
  dipoleTransform.translate(gOpt.dipoleTranslation);

  auto dipoleNode = std::make_unique<Acts::Experimental::Blueprint::Node>(
      "Dipole", dipoleTransform, Acts::VolumeBounds::eCuboid, gOpt.dipoleBounds,
      dipoleLayerBuilder);

  trackerBP->add(std::move(dipoleNode));

  // Quadrupole 1
  auto quad1ZBounds =
      std::make_tuple(gOpt.quad1Translation.z() - gOpt.quad1Bounds.at(2),
                      gOpt.quad1Translation.z() + gOpt.quad1Bounds.at(2));

  auto quad1LayerBuilder = makeLayerBuilder(
      world, gOpt.g4ToWorld, names, {quad1ZBounds}, {Acts::BinningValue::binZ});

  Acts::Transform3 quad1Transform = Acts::Transform3::Identity();
  quad1Transform.rotate(gOpt.actsToWorld.rotation().inverse());
  quad1Transform.translate(gOpt.quad1Translation);

  auto quad1Node = std::make_unique<Acts::Experimental::Blueprint::Node>(
      "Quad1", quad1Transform, Acts::VolumeBounds::eCuboid, gOpt.quad1Bounds,
      quad1LayerBuilder);

  trackerBP->add(std::move(quad1Node));

  // Quadrupole 2
  auto quad2ZBounds =
      std::make_tuple(gOpt.quad2Translation.z() - gOpt.quad2Bounds.at(2),
                      gOpt.quad2Translation.z() + gOpt.quad2Bounds.at(2));

  auto quad2LayerBuilder = makeLayerBuilder(
      world, gOpt.g4ToWorld, names, {quad2ZBounds}, {Acts::BinningValue::binZ});

  Acts::Transform3 quad2Transform = Acts::Transform3::Identity();
  quad2Transform.rotate(gOpt.actsToWorld.rotation().inverse());
  quad2Transform.translate(gOpt.quad2Translation);

  auto quad2Node = std::make_unique<Acts::Experimental::Blueprint::Node>(
      "Quad2", quad2Transform, Acts::VolumeBounds::eCuboid, gOpt.quad2Bounds,
      quad2LayerBuilder);

  trackerBP->add(std::move(quad2Node));

  // Quadrupole 3
  auto quad3ZBounds =
      std::make_tuple(gOpt.quad3Translation.z() - gOpt.quad3Bounds.at(2),
                      gOpt.quad3Translation.z() + gOpt.quad3Bounds.at(2));

  auto quad3LayerBuilder = makeLayerBuilder(
      world, gOpt.g4ToWorld, names, {quad3ZBounds}, {Acts::BinningValue::binZ});

  Acts::Transform3 quad3Transform = Acts::Transform3::Identity();
  quad3Transform.rotate(gOpt.actsToWorld.rotation().inverse());
  quad3Transform.translate(gOpt.quad3Translation);

  auto quad3Node = std::make_unique<Acts::Experimental::Blueprint::Node>(
      "Quad3", quad3Transform, Acts::VolumeBounds::eCuboid, gOpt.quad3Bounds,
      quad3LayerBuilder);

  trackerBP->add(std::move(quad3Node));

  return trackerBP;
};

std::shared_ptr<const Acts::Experimental::Detector> buildE320Detector(
    const std::unique_ptr<Acts::Experimental::Blueprint::Node> detectorBpr,
    const Acts::GeometryContext& gctx,
    const E320Geometry::GeometryOptions& gOpt,
    const std::vector<Acts::GeometryIdentifier>& materialVetos) {
  // Complete and fill gaps
  Acts::Experimental::detail::BlueprintHelper::fillGaps(*detectorBpr, false);

  auto detectorBuilder =
      std::make_shared<Acts::Experimental::CuboidalContainerBuilder>(
          *detectorBpr, Acts::Logging::VERBOSE);

  auto idGenCfg = E320GeometryIdGenerator::Config{false, 0u, true, false, gOpt};

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
  dCfg.geoIdGenerator = std::make_shared<E320GeometryIdGenerator>(
      idGenCfg,
      Acts::getDefaultLogger("GeoIdGenerator", Acts::Logging::VERBOSE));

  auto detector = Acts::Experimental::DetectorBuilder(dCfg).construct(gctx);

  return detector;
}

std::shared_ptr<const Acts::Experimental::Detector> buildE320Detector(
    const std::unique_ptr<Acts::Experimental::Blueprint::Node> detectorBpr,
    const Acts::GeometryContext& gctx,
    const E320Geometry::GeometryOptions& gOpt,
    const std::string jsonMaterialPath,
    const std::vector<Acts::GeometryIdentifier>& materialVetos) {
  // Complete and fill gaps
  Acts::Experimental::detail::BlueprintHelper::fillGaps(*detectorBpr, false);

  auto detectorBuilder =
      std::make_shared<Acts::Experimental::CuboidalContainerBuilder>(
          *detectorBpr, Acts::Logging::VERBOSE);

  auto idGenCfg = E320GeometryIdGenerator::Config{false, 0u, true, false, gOpt};

  // Material provider
  auto jMatMapConverterCfg = Acts::MaterialMapJsonConverter::Config();
  jMatMapConverterCfg.context = gctx;

  auto jMatDecorator = std::make_shared<Acts::JsonMaterialDecorator>(
      jMatMapConverterCfg, jsonMaterialPath, Acts::Logging::VERBOSE,
      materialVetos);

  // Detector builder
  Acts::Experimental::DetectorBuilder::Config dCfg;
  dCfg.auxiliary = "Detector builder";
  dCfg.name = "Detector from blueprint";
  dCfg.builder = detectorBuilder;
  dCfg.materialDecorator = jMatDecorator;
  dCfg.geoIdGenerator = std::make_shared<E320GeometryIdGenerator>(
      idGenCfg,
      Acts::getDefaultLogger("GeoIdGenerator", Acts::Logging::VERBOSE));

  auto detector = Acts::Experimental::DetectorBuilder(
                      dCfg, Acts::getDefaultLogger("DetectorBuilder",
                                                   Acts::Logging::VERBOSE))
                      .construct(gctx);
  return detector;
}

}  // namespace E320Geometry
