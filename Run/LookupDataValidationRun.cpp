#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <memory>

#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/DummyReader.hpp"
#include "TrackingPipeline/Io/RootTrackLookupValidationWriter.hpp"
#include "TrackingPipeline/MagneticField/CompositeMagField.hpp"
#include "TrackingPipeline/MagneticField/ConstantBoundedField.hpp"
#include "TrackingPipeline/MagneticField/DipoleMagField.hpp"
#include "TrackingPipeline/MagneticField/QuadrupoleMagField.hpp"
#include "TrackingPipeline/Simulation/BremsstrahlungMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/GaussianVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"
#include "TrackingPipeline/Simulation/SimpleDigitizer.hpp"
#include "TrackingPipeline/TrackFinding/DipoleTrackLookupProvider.hpp"
#include "TrackingPipeline/TrackFinding/TrackLookupValidationAlgorithm.hpp"

using namespace Acts::UnitLiterals;

using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                    Acts::Experimental::DetectorNavigator>;
using TrackParameters = Acts::CurvilinearTrackParameters;

/// @brief Run the propagation through
/// a uniform energy spectrum and record the
/// energy vs position histograms for each layer
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
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_gdmls/"
      "ett_geometry_f566a577.gdml";
  std::vector<std::string> names{"OPPPSensitive", "DetChamberWindow"};

  // Veto PDC window material mapping
  // to preserve homogeneous material
  // from Geant4
  Acts::GeometryIdentifier pdcWindowId;
  pdcWindowId.setApproach(1);
  std::vector<Acts::GeometryIdentifier> materialVeto{pdcWindowId};

  std::string materialPath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Prototype/"
      "E320Prototype_material/"
      "Uniform_DirectZ_TrackerOnly_256x128_1M/material.json";

  // Build the detector
  auto trackerBP = E320Geometry::makeBlueprintE320(gdmlPath, names, gOpt);
  auto detector = E320Geometry::buildE320Detector(
      std::move(trackerBP), gctx, gOpt, materialPath, materialVeto);

  // --------------------------------------------------------------
  // The magnetic field setup

  // Extent in already rotated frame
  Acts::Extent quad1Extent;
  quad1Extent.set(Acts::BinningValue::binX,
                  gOpt.quad1Translation[0] - gOpt.quad1Bounds[0],
                  gOpt.quad1Translation[0] + gOpt.quad1Bounds[0]);
  quad1Extent.set(Acts::BinningValue::binZ,
                  gOpt.quad1Translation[1] - gOpt.quad1Bounds[1],
                  gOpt.quad1Translation[1] + gOpt.quad1Bounds[1]);
  quad1Extent.set(Acts::BinningValue::binY,
                  gOpt.quad1Translation[2] - gOpt.quad1Bounds[2],
                  gOpt.quad1Translation[2] + gOpt.quad1Bounds[2]);

  Acts::Extent quad2Extent;
  quad2Extent.set(Acts::BinningValue::binX,
                  gOpt.quad2Translation[0] - gOpt.quad2Bounds[0],
                  gOpt.quad2Translation[0] + gOpt.quad2Bounds[0]);
  quad2Extent.set(Acts::BinningValue::binZ,
                  gOpt.quad2Translation[1] - gOpt.quad2Bounds[1],
                  gOpt.quad2Translation[1] + gOpt.quad2Bounds[1]);
  quad2Extent.set(Acts::BinningValue::binY,
                  gOpt.quad2Translation[2] - gOpt.quad2Bounds[2],
                  gOpt.quad2Translation[2] + gOpt.quad2Bounds[2]);

  Acts::Extent quad3Extent;
  quad3Extent.set(Acts::BinningValue::binX,
                  gOpt.quad3Translation[0] - gOpt.quad3Bounds[0],
                  gOpt.quad3Translation[0] + gOpt.quad3Bounds[0]);
  quad3Extent.set(Acts::BinningValue::binZ,
                  gOpt.quad3Translation[1] - gOpt.quad3Bounds[1],
                  gOpt.quad3Translation[1] + gOpt.quad3Bounds[1]);
  quad3Extent.set(Acts::BinningValue::binY,
                  gOpt.quad3Translation[2] - gOpt.quad3Bounds[2],
                  gOpt.quad3Translation[2] + gOpt.quad3Bounds[2]);

  Acts::Extent dipoleExtent;
  dipoleExtent.set(Acts::BinningValue::binX,
                   gOpt.dipoleTranslation.x() - gOpt.dipoleBounds[0],
                   gOpt.dipoleTranslation.x() + gOpt.dipoleBounds[0]);
  dipoleExtent.set(Acts::BinningValue::binZ,
                   gOpt.dipoleTranslation.y() - gOpt.dipoleBounds[1],
                   gOpt.dipoleTranslation.y() + gOpt.dipoleBounds[1]);
  dipoleExtent.set(Acts::BinningValue::binY,
                   gOpt.dipoleTranslation.z() - gOpt.dipoleBounds[2],
                   gOpt.dipoleTranslation.z() + gOpt.dipoleBounds[2]);

  Acts::Extent xCorrectorExtent;
  xCorrectorExtent.set(
      Acts::BinningValue::binX,
      gOpt.xCorrectorTranslation.x() - gOpt.xCorrectorBounds[0],
      gOpt.xCorrectorTranslation.x() + gOpt.xCorrectorBounds[0]);
  xCorrectorExtent.set(
      Acts::BinningValue::binZ,
      gOpt.xCorrectorTranslation.y() - gOpt.xCorrectorBounds[1],
      gOpt.xCorrectorTranslation.y() + gOpt.xCorrectorBounds[1]);
  xCorrectorExtent.set(
      Acts::BinningValue::binY,
      gOpt.xCorrectorTranslation.z() - gOpt.xCorrectorBounds[2],
      gOpt.xCorrectorTranslation.z() + gOpt.xCorrectorBounds[2]);

  QuadrupoleMagField quad1Field(
      gOpt.quadrupolesParams[0],
      gOpt.actsToWorldRotation.inverse() * gOpt.quad1Translation,
      gOpt.actsToWorldRotation);
  QuadrupoleMagField quad2Field(
      gOpt.quadrupolesParams[1],
      gOpt.actsToWorldRotation.inverse() * gOpt.quad2Translation,
      gOpt.actsToWorldRotation);
  QuadrupoleMagField quad3Field(
      gOpt.quadrupolesParams[2],
      gOpt.actsToWorldRotation.inverse() * gOpt.quad3Translation,
      gOpt.actsToWorldRotation);

  double dipoleB = 0.2192_T;
  DipoleMagField dipoleField(
      gOpt.dipoleParams, dipoleB, gOpt.actsToWorldRotation,
      gOpt.actsToWorldRotation.inverse() * gOpt.dipoleTranslation);

  Acts::Vector3 xCorrectorB(0, 0, -0.026107_T);
  ConstantBoundedField xCorrectorField(xCorrectorB, xCorrectorExtent);

  CompositeMagField::FieldComponents fieldComponents = {
      {quad1Extent, &quad1Field},
      {quad2Extent, &quad2Field},
      {quad3Extent, &quad3Field},
      {dipoleExtent, &dipoleField},
      {xCorrectorExtent, &xCorrectorField}};

  auto field = std::make_shared<CompositeMagField>(fieldComponents);

  auto aStore =
      std::make_shared<std::map<Acts::GeometryIdentifier, Acts::Transform3>>();

  double detectorTilt = -0.00297;
  std::map<int, Acts::Vector3> shifts{
      {8, Acts::Vector3(-11.1_mm, -0.75_mm, 0_mm)},
      {6, Acts::Vector3(-11.1_mm, -0.75_mm, 0_mm)},
      {4, Acts::Vector3(-11.1_mm, -0.75_mm, 0_mm)},
      {2, Acts::Vector3(-11.1_mm, -0.75_mm, 0_mm)},
      {0, Acts::Vector3(-11.1_mm, -0.75_mm, 0_mm)}};

  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive()) {
        Acts::Transform3 nominal = s->transform(Acts::GeometryContext());
        nominal.pretranslate(shifts.at(s->geometryId().sensitive() - 1));
        aStore->emplace(s->geometryId(), nominal);
      }
    }
  }
  AlignmentContext alignCtx(aStore);
  gctx = Acts::GeometryContext{alignCtx};

  Acts::GeometryIdentifier midGeoId;
  midGeoId.setSensitive(5);
  Acts::Vector3 detectorCenter = detector->findSurface(midGeoId)->center(gctx);
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive()) {
        Acts::Transform3 nominal = aStore->at(s->geometryId());
        nominal.pretranslate(-detectorCenter);
        nominal.prerotate(Acts::AngleAxis3(detectorTilt, Acts::Vector3::UnitX())
                              .toRotationMatrix());
        nominal.pretranslate(detectorCenter);
        aStore->at(s->geometryId()) = nominal;
      }
    }
  }

  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive()) {
        std::cout << "-----------------------------------\n";
        std::cout << "SURFACE " << s->geometryId() << "\n";
        std::cout << "CENTER " << s->center(gctx).transpose() << " -- "
                  << s->center(Acts::GeometryContext()).transpose() << "\n";
        std::cout << "NORMAL "
                  << s->normal(gctx, s->center(gctx), Acts::Vector3::UnitY())
                         .transpose()
                  << " -- "
                  << s->normal(Acts::GeometryContext(),
                               s->center(Acts::GeometryContext()),
                               Acts::Vector3::UnitY())
                         .transpose()
                  << "\n";
        std::cout << "ROTATION \n"
                  << s->transform(gctx).rotation() << " -- \n"
                  << "\n"
                  << s->transform(Acts::GeometryContext()).rotation() << "\n";
        std::cout << "EXTENT "
                  << s->polyhedronRepresentation(gctx, 1000).extent()
                  << "\n -- \n"
                  << s->polyhedronRepresentation(Acts::GeometryContext(), 1000)
                         .extent()
                  << "\n";
      }
    }
  }

  // --------------------------------------------------------------
  // Event reading
  SimpleSourceLink::SurfaceAccessor surfaceAccessor{detector.get()};

  // Setup the sequencer
  Sequencer::Config seqCfg;
  seqCfg.numThreads = 1;
  seqCfg.trackFpes = false;
  seqCfg.logLevel = logLevel;
  Sequencer sequencer(seqCfg);

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));

  // --------------------------------------------------------------
  // Add dummy reader
  DummyReader::Config dummyReaderCfg;
  dummyReaderCfg.outputSourceLinks = "SimMeasurements";
  dummyReaderCfg.outputSimClusters = "SimClusters";
  dummyReaderCfg.nEvents = 1e4;

  sequencer.addReader(std::make_shared<DummyReader>(dummyReaderCfg));

  // --------------------------------------------------------------
  // Simulate track propagation

  // Setup the measurements creator
  Acts::Experimental::DetectorNavigator::Config cptNavCfg;
  cptNavCfg.detector = detector.get();
  cptNavCfg.resolvePassive = false;
  cptNavCfg.resolveMaterial = true;
  cptNavCfg.resolveSensitive = true;

  Acts::Experimental::DetectorNavigator measCreatorNavigator(
      cptNavCfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));
  Acts::EigenStepper<> measCreatorStepper(field);

  Propagator measCreatorPropagator(std::move(measCreatorStepper),
                                   std::move(measCreatorNavigator));

  // Digitizer
  SimpleDigitizer::Config digitizerCfg;
  digitizerCfg.resolution = {5_um, 5_um};
  auto digitizer = std::make_shared<SimpleDigitizer>(digitizerCfg);

  // Vertex generator
  auto vertexGen = std::make_shared<GaussianVertexGenerator>(
      Acts::Vector3(0, gOpt.beWindowTranslation[2], 0),
      75_um * 75_um * Acts::SquareMatrix3::Identity());

  // Momentum generator
  auto momGen = std::make_shared<BremsstrahlungMomentumGenerator>(
      Acts::Vector3(0, 0, 0),
      0.3_MeV * 0.3_MeV * Acts::SquareMatrix3::Identity());

  // Measurement creator
  MeasurementsCreator::Config measCreatorCfg;
  measCreatorCfg.vertexGenerator = vertexGen;
  measCreatorCfg.momentumGenerator = momGen;
  measCreatorCfg.hitDigitizer = digitizer;
  measCreatorCfg.maxSteps = 1000;
  measCreatorCfg.isSignal = true;

  auto measCreator = std::make_shared<MeasurementsCreator>(
      measCreatorPropagator, measCreatorCfg);

  MeasurementsEmbeddingAlgorithm::Config measCreatorAlgoCfg;
  measCreatorAlgoCfg.inputSourceLinks = "SimMeasurements";
  measCreatorAlgoCfg.inputSimClusters = "SimClusters";
  measCreatorAlgoCfg.outputSourceLinks = "Measurements";
  measCreatorAlgoCfg.outputSimClusters = "Clusters";
  measCreatorAlgoCfg.measurementGenerator = measCreator;
  measCreatorAlgoCfg.randomNumberSvc =
      std::make_shared<RandomNumbers>(RandomNumbers::Config());
  measCreatorAlgoCfg.nMeasurements = 1;

  sequencer.addAlgorithm(std::make_shared<MeasurementsEmbeddingAlgorithm>(
      measCreatorAlgoCfg, logLevel));

  // --------------------------------------------------------------
  // Lookup data validation

  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> refLayers;
  const auto& refVolume = detector->findDetectorVolume("layer0");
  for (const auto* surf : refVolume->surfaces()) {
    refLayers.try_emplace(surf->geometryId(), surf);
  }

  // Validation algorithm
  TrackLookupValidationAlgorithm::Config validatorCfg;

  E320DipoleTrackLookupProvider::Config lookupProviderCfg;
  lookupProviderCfg.dipoleAmplidute = 0.2192;
  lookupProviderCfg.dipolePosition = gOpt.dipoleTranslation[2];
  lookupProviderCfg.dipoleSize = 0.914;

  lookupProviderCfg.correctorAmplidute = -0.026107_T;
  lookupProviderCfg.correctorPosition = gOpt.xCorrectorTranslation[2];
  lookupProviderCfg.correctorSize = 0.23368;

  lookupProviderCfg.layerPosition = gOpt.staveZ.at(8);
  lookupProviderCfg.referenceSurface = refLayers.begin()->second;
  E320DipoleTrackLookupProvider lookupProvider(lookupProviderCfg);

  validatorCfg.refLayers = refLayers;
  validatorCfg.estimator.connect<&E320DipoleTrackLookupProvider::lookup>(
      &lookupProvider);
  validatorCfg.inputClusters = "Clusters";
  validatorCfg.outputIpPars = "ipPars";
  validatorCfg.outputIpParsEst = "ipParsEst";
  validatorCfg.outputRefLayerPars = "refPars";
  validatorCfg.outputRefLayerParsEst = "refParsEst";

  sequencer.addAlgorithm(
      std::make_shared<TrackLookupValidationAlgorithm>(validatorCfg, logLevel));

  RootTrackLookupValidationWriter::Config validationWriterCfg;
  validationWriterCfg.inputIpPars = "ipPars";
  validationWriterCfg.inputIpParsEst = "ipParsEst";
  validationWriterCfg.inputRefLayerPars = "refPars";
  validationWriterCfg.inputRefLayerParsEst = "refParsEst";
  validationWriterCfg.path = "lookup-validation.root";
  validationWriterCfg.treeName = "validation";

  sequencer.addWriter(std::make_shared<RootTrackLookupValidationWriter>(
      validationWriterCfg, logLevel));

  return sequencer.run();
}
