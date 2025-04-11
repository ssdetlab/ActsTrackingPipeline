#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <filesystem>
#include <initializer_list>
#include <memory>

#include "TrackingPipeline/Clustering/HourglassFilter.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/DummyReader.hpp"
#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"
#include "TrackingPipeline/Io/RootTrackParamsReader.hpp"
#include "TrackingPipeline/MagneticField/CompositeMagField.hpp"
#include "TrackingPipeline/MagneticField/DipoleMagField.hpp"
#include "TrackingPipeline/MagneticField/QuadrupoleMagField.hpp"
#include "TrackingPipeline/Simulation/E320BeamBkgVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/E320CptBkgGenerator.hpp"
#include "TrackingPipeline/Simulation/E320CptBkgVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/E320HistDigitizer.hpp"
#include "TrackingPipeline/Simulation/KDEMomentumxGenerator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"

using ActionList = Acts::ActionList<>;
using AbortList = Acts::AbortList<Acts::EndOfWorldReached>;

using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                    Acts::Experimental::DetectorNavigator>;
using PropagatorOptions =
    typename Propagator::template Options<ActionList, AbortList>;

using KFTrajectory = Acts::VectorMultiTrajectory;
using KFTrackContainer = Acts::VectorTrackContainer;
using KF = Acts::KalmanFitter<Propagator, KFTrajectory>;

using CKFTrackContainer = Acts::TrackContainer<Acts::VectorTrackContainer,
                                               Acts::VectorMultiTrajectory,
                                               Acts::detail::ValueHolder>;

using TrackStateContainerBackend =
    typename CKFTrackContainer::TrackStateContainerBackend;

using namespace Acts::UnitLiterals;

// TODO: There's something wrong with the KDE
int main() {
  // Set the log level
  Acts::Logging::Level logLevel = Acts::Logging::FATAL;

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
      "full_detector_sim/ettgeom_magnet_pdc_tracker.gdml";
  std::vector<std::string> names{"OPPPSensitive", "DetChamberWindow"};

  // Veto PDC window material mapping
  // to preserve homogeneous material
  // from Geant4
  Acts::GeometryIdentifier pdcWindowId;
  pdcWindowId.setApproach(1);
  std::vector<Acts::GeometryIdentifier> materialVeto{pdcWindowId};

  std::string materialPath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_material/"
      "full_detector_sim/Uniform_DirectZ_TrackerOnly_256x128_1M/material.json";

  // Build the detector
  auto trackerBP = E320Geometry::makeBlueprintE320(gdmlPath, names, gOpt);
  auto detector = E320Geometry::buildE320Detector(
      std::move(trackerBP), gctx, gOpt, materialPath, materialVeto);

  for (auto& v : detector->volumes()) {
    std::cout << v->name() << std::endl;
    for (auto& s : v->surfaces()) {
      std::cout << s->center(gctx).transpose() << "   "
                << s->normal(gctx, s->center(gctx), Acts::Vector3(1, 0, 0))
                       .transpose()
                << std::endl;
    }
  }

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

  double dipoleB = 0.31_T;
  DipoleMagField dipoleField(
      gOpt.dipoleParams, dipoleB, gOpt.actsToWorldRotation,
      gOpt.actsToWorldRotation.inverse() * gOpt.dipoleTranslation);

  CompositeMagField::FieldComponents fieldComponents = {
      {quad1Extent, &quad1Field},
      {quad2Extent, &quad2Field},
      {quad3Extent, &quad3Field},
      {dipoleExtent, &dipoleField}};

  auto field = std::make_shared<CompositeMagField>(fieldComponents);

  // --------------------------------------------------------------
  // Event reading

  SimpleSourceLink::SurfaceAccessor surfaceAccessor{detector.get()};
  auto hourglassFilter = std::make_shared<HourglassFilter>();
  hourglassFilter->surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  // Setup the sequencer
  Sequencer::Config seqCfg;
  seqCfg.numThreads = 1;
  seqCfg.trackFpes = false;
  Sequencer sequencer(seqCfg);

  // --------------------------------------------------------------
  // Add dummy reader
  DummyReader::Config dummyReaderCfg;
  dummyReaderCfg.outputSourceLinks = "SimMeasurements";
  dummyReaderCfg.outputSimClusters = "SimClusters";
  dummyReaderCfg.nEvents = 1e1;

  sequencer.addReader(std::make_shared<DummyReader>(dummyReaderCfg));

  // --------------------------------------------------------------
  // Compton background embedding

  // Setup the measurements creator
  Acts::Experimental::DetectorNavigator::Config cptNavCfg;
  cptNavCfg.detector = detector.get();
  cptNavCfg.resolvePassive = false;
  cptNavCfg.resolveMaterial = true;
  cptNavCfg.resolveSensitive = true;

  Acts::Experimental::DetectorNavigator cptNavigator(
      cptNavCfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));
  Acts::EigenStepper<> cptStepper(field);

  Propagator cptBkgPropagator(std::move(cptStepper), std::move(cptNavigator));

  // Digitizer
  E320Sim::E320HistDigitizer::Config cptDigitizerCfg;
  cptDigitizerCfg.pathToHist =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_genHists/"
      "comptonBkg/genHistSize.root";
  cptDigitizerCfg.histName = "histSize";

  auto cptDigitizer =
      std::make_shared<E320Sim::E320HistDigitizer>(cptDigitizerCfg);

  // Track parameters generator
  RootTrackParamsReader::Config cptReaderCfg;
  cptReaderCfg.treeName = "track-parameters";
  cptReaderCfg.transform = Acts::Transform3::Identity();
  std::string pathToDirCpt =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_genRootData/"
      "comptonBkg";

  // Get the paths to the files in the directory
  for (const auto& entry : std::filesystem::directory_iterator(pathToDirCpt)) {
    std::string pathToFile = entry.path();
    cptReaderCfg.filePaths.push_back(pathToFile);
  }

  E320Sim::E320CptBkgGenerator::Config cptBkgGenCfg;
  cptBkgGenCfg.nIterations = 1;
  cptBkgGenCfg.sensitivity = 0.5;
  cptBkgGenCfg.vertexGenCfg = E320Sim::cptBkgVertexGenConfig();
  cptBkgGenCfg.trackParamsReader =
      std::make_shared<RootTrackParamsReader>(cptReaderCfg);
  auto cptBkgGen = std::make_shared<E320Sim::E320CptBkgGenerator>(cptBkgGenCfg);

  // Measurement creator
  MeasurementsCreator::Config cptMcCfg;
  cptMcCfg.vertexGenerator = cptBkgGen;
  cptMcCfg.momentumGenerator = cptBkgGen;
  cptMcCfg.hitDigitizer = cptDigitizer;
  cptMcCfg.maxSteps = 100;
  cptMcCfg.isSignal = false;

  auto cptMeasurementsCreator =
      std::make_shared<MeasurementsCreator>(cptBkgPropagator, cptMcCfg);

  MeasurementsEmbeddingAlgorithm::Config cptMceCfg;
  cptMceCfg.inputSourceLinks = "SimMeasurements";
  cptMceCfg.inputSimClusters = "SimClusters";
  cptMceCfg.outputSourceLinks = "Measurements";
  cptMceCfg.outputSimClusters = "Clusters";
  cptMceCfg.measurementGenerator = cptMeasurementsCreator;
  //   cptMceCfg.clusterFilter = hourglassFilter;
  cptMceCfg.randomNumberSvc =
      std::make_shared<RandomNumbers>(RandomNumbers::Config());
  cptMceCfg.nMeasurements = 660;

  sequencer.addAlgorithm(
      std::make_shared<MeasurementsEmbeddingAlgorithm>(cptMceCfg, logLevel));

  // --------------------------------------------------------------
  // Beam background embedding

  // Setup the measurements creator
  Acts::Experimental::DetectorNavigator::Config beamNavCfg;
  beamNavCfg.detector = detector.get();
  beamNavCfg.resolvePassive = false;
  beamNavCfg.resolveMaterial = true;
  beamNavCfg.resolveSensitive = true;

  Acts::Experimental::DetectorNavigator beamNavigator(
      beamNavCfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));
  Acts::EigenStepper<> beamStepper(field);

  Propagator beamBkgPropagator(std::move(beamStepper),
                               std::move(beamNavigator));

  // Digitizer
  E320Sim::E320HistDigitizer::Config beamDigitizerCfg;
  beamDigitizerCfg.pathToHist =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_genHists/beamBkg/"
      "genHistSize.root";
  beamDigitizerCfg.histName = "histSize";

  auto beamDigitizer =
      std::make_shared<E320Sim::E320HistDigitizer>(beamDigitizerCfg);

  // Vertex generator
  auto beamVertexGen = std::make_shared<E320Sim::E320BeamBkgVertexGenerator>();

  beamVertexGen->yBoundLow = gOpt.chipY.at(0) - gOpt.chipSizeY / 2;
  beamVertexGen->yBoundHigh = gOpt.chipY.at(8) + gOpt.chipSizeY / 2;

  beamVertexGen->xBoundLow = gOpt.chipX - gOpt.chipSizeX / 2;
  beamVertexGen->xBoundHigh = gOpt.chipX + gOpt.chipSizeX / 2;

  beamVertexGen->yPower = 69.8048;
  beamVertexGen->yShift = -22165;
  beamVertexGen->zProbs = {0.25, 0.24, 0.23, 0.28};
  beamVertexGen->zPositions = {gOpt.staveZ.at(0), gOpt.staveZ.at(1),
                               gOpt.staveZ.at(2), gOpt.staveZ.at(3)};

  // Momentum generator
  RootTrackParamsReader::Config beamReaderCfg;
  beamReaderCfg.treeName = "track-parameters";
  beamReaderCfg.transform = Acts::Transform3::Identity();
  std::string pathToDirBeam =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_genRootData/"
      "beamBkg";

  // Get the paths to the files in the directory
  for (const auto& entry : std::filesystem::directory_iterator(pathToDirBeam)) {
    std::string pathToFile = entry.path();
    beamReaderCfg.filePaths.push_back(pathToFile);
  }

  KDEMomentumGenerator::Config beamMomGenCfg;
  beamMomGenCfg.trackParamsReader =
      std::make_shared<RootTrackParamsReader>(beamReaderCfg);
  beamMomGenCfg.nIterations = 1;
  beamMomGenCfg.sensitivity = 0.5;
  beamMomGenCfg.transform =
      Acts::Transform3(Acts::Translation3(Acts::Vector3(0, 0, 0)) *
                       gOpt.actsToWorldRotation.inverse());

  auto beamMomGen = std::make_shared<KDEMomentumGenerator>(beamMomGenCfg);

  // Measurement creator
  MeasurementsCreator::Config beamMcCfg;
  beamMcCfg.vertexGenerator = beamVertexGen;
  beamMcCfg.momentumGenerator = beamMomGen;
  beamMcCfg.hitDigitizer = beamDigitizer;
  beamMcCfg.maxSteps = 100;
  beamMcCfg.isSignal = false;

  auto beamMeasurementsCreator =
      std::make_shared<MeasurementsCreator>(beamBkgPropagator, beamMcCfg);

  MeasurementsEmbeddingAlgorithm::Config beamMceCfg;
  beamMceCfg.inputSourceLinks = "SimMeasurements";
  beamMceCfg.inputSimClusters = "SimClusters";
  beamMceCfg.outputSourceLinks = "Measurements";
  beamMceCfg.outputSimClusters = "Clusters";
  beamMceCfg.measurementGenerator = beamMeasurementsCreator;
  //   beamMceCfg.clusterFilter = hourglassFilter;
  beamMceCfg.randomNumberSvc =
      std::make_shared<RandomNumbers>(RandomNumbers::Config());
  beamMceCfg.nMeasurements = 3450;

  sequencer.addAlgorithm(
      std::make_shared<MeasurementsEmbeddingAlgorithm>(beamMceCfg, logLevel));

  // --------------------------------------------------------------
  // Event write out

  // Sim cluster writer
  auto clusterWriterCfg = RootSimClusterWriter::Config();

  clusterWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  clusterWriterCfg.inputClusters = "Clusters";
  clusterWriterCfg.treeName = "clusters";
  clusterWriterCfg.filePath = "clusters-ncs-beam.root";

  sequencer.addWriter(
      std::make_shared<RootSimClusterWriter>(clusterWriterCfg, logLevel));

  return sequencer.run();
}
