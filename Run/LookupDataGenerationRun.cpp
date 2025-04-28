#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/DummyReader.hpp"
#include "TrackingPipeline/Io/JsonTrackLookupWriter.hpp"
#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"
#include "TrackingPipeline/MagneticField/CompositeMagField.hpp"
#include "TrackingPipeline/MagneticField/ConstantBoundedField.hpp"
#include "TrackingPipeline/MagneticField/DipoleMagField.hpp"
#include "TrackingPipeline/MagneticField/QuadrupoleMagField.hpp"
#include "TrackingPipeline/Simulation/IdealDigitizer.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"
#include "TrackingPipeline/Simulation/RangedUniformMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/StationaryVertexGenerator.hpp"
#include "TrackingPipeline/TrackFinding/TrackLookupEstimationAlgorithm.hpp"

using namespace Acts::UnitLiterals;

using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                    Acts::Experimental::DetectorNavigator>;
using TrackParameters = Acts::CurvilinearTrackParameters;

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

  // TODO: Check if it's the real field
  // Acts::Vector3 xCorrectorB(0, 0, -0.0536_T);
  Acts::Vector3 xCorrectorB(0, 0, 0);
  ConstantBoundedField xCorrectorField(xCorrectorB, xCorrectorExtent);

  CompositeMagField::FieldComponents fieldComponents = {
      {quad1Extent, &quad1Field},
      {quad2Extent, &quad2Field},
      {quad3Extent, &quad3Field},
      {xCorrectorExtent, &xCorrectorField},
      {dipoleExtent, &dipoleField}};

  auto field = std::make_shared<CompositeMagField>(fieldComponents);

  // --------------------------------------------------------------
  // Event reading
  SimpleSourceLink::SurfaceAccessor surfaceAccessor{detector.get()};

  // Setup the sequencer
  Sequencer::Config seqCfg;
  // seqCfg.events = 1;
  seqCfg.numThreads = 1;
  seqCfg.trackFpes = false;
  Sequencer sequencer(seqCfg);

  // --------------------------------------------------------------
  // Add dummy reader
  DummyReader::Config dummyReaderCfg;
  dummyReaderCfg.outputSourceLinks = "SimMeasurements";
  dummyReaderCfg.outputSimClusters = "SimClusters";
  dummyReaderCfg.nEvents = 2e5;

  sequencer.addReader(std::make_shared<DummyReader>(dummyReaderCfg));

  // --------------------------------------------------------------
  // Setup the measurements creator

  Acts::Experimental::DetectorNavigator::Config navCfg;
  navCfg.detector = detector.get();
  navCfg.resolvePassive = false;
  navCfg.resolveMaterial = true;
  navCfg.resolveSensitive = true;

  Acts::Experimental::DetectorNavigator navigator(
      navCfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));
  Acts::EigenStepper<> stepper(field);

  auto propagator = Propagator(std::move(stepper), std::move(navigator));

  auto momGen = std::make_shared<RangedUniformMomentumGenerator>();
  momGen->Pranges = {{1.9_GeV, 2.1_GeV},
                     {2.1_GeV, 2.3_GeV},
                     {2.3_GeV, 2.5_GeV},
                     {2.5_GeV, 2.7_GeV},
                     {2.7_GeV, 2.9_GeV}};

  auto vertexGen = std::make_shared<StationaryVertexGenerator>();
  vertexGen->vertex =
      Acts::Vector3(gOpt.beWindowTranslation[0], gOpt.beWindowTranslation[2],
                    -gOpt.beWindowTranslation[1]);

  // Digitizer
  auto digitizer = std::make_shared<IdealDigitizer>();

  MeasurementsCreator::Config mcCfg;
  mcCfg.vertexGenerator = vertexGen;
  mcCfg.momentumGenerator = momGen;
  mcCfg.hitDigitizer = digitizer;
  mcCfg.maxSteps = 300;

  auto measurementsCreator =
      std::make_shared<MeasurementsCreator>(propagator, mcCfg);

  MeasurementsEmbeddingAlgorithm::Config mcaCfg;
  mcaCfg.inputSourceLinks = "SimMeasurements";
  mcaCfg.inputSimClusters = "SimClusters";
  mcaCfg.outputSourceLinks = "Measurements";
  mcaCfg.outputSimClusters = "Clusters";
  mcaCfg.measurementGenerator = measurementsCreator;
  mcaCfg.randomNumberSvc =
      std::make_shared<RandomNumbers>(RandomNumbers::Config());
  mcaCfg.nMeasurements = 1;

  sequencer.addAlgorithm(
      std::make_shared<MeasurementsEmbeddingAlgorithm>(mcaCfg, logLevel));

  // --------------------------------------------------------------
  // Lookup data generation

  JsonTrackLookupWriter::Config lookupWriterCfg;
  lookupWriterCfg.path = "lookup-prototype.json";

  auto lookupWriter = std::make_shared<JsonTrackLookupWriter>(lookupWriterCfg);

  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> refLayers;
  const auto& refVolume = detector->findDetectorVolume("layer0");
  for (const auto* surf : refVolume->surfaces()) {
    refLayers.try_emplace(surf->geometryId(), surf);
  }

  TrackLookupEstimationAlgorithm::Config estimatorCfg;
  estimatorCfg.trackLookupGridWriters = {lookupWriter};
  estimatorCfg.refLayers = refLayers;
  estimatorCfg.bins = {1000, 1};
  estimatorCfg.inputClusters = "Clusters";

  sequencer.addAlgorithm(
      std::make_shared<TrackLookupEstimationAlgorithm>(estimatorCfg, logLevel));

  // --------------------------------------------------------------
  // Event write out
  //
  // Sim cluster writer
  auto clusterWriterCfg = RootSimClusterWriter::Config();

  clusterWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  clusterWriterCfg.inputClusters = "Clusters";
  clusterWriterCfg.treeName = "clusters";
  clusterWriterCfg.filePath = "clusters-protoype.root";

  sequencer.addWriter(
      std::make_shared<RootSimClusterWriter>(clusterWriterCfg, logLevel));
  return sequencer.run();
}
