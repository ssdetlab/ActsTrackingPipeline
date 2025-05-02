#include "Acts/Utilities/Logger.hpp"

#include <filesystem>
#include <memory>

#include "TrackingPipeline/Clustering/HourglassFilter.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/DummyReader.hpp"
#include "TrackingPipeline/Io/E320RootSimDataReader.hpp"
#include "TrackingPipeline/Io/JsonTrackLookupReader.hpp"
#include "TrackingPipeline/Io/JsonTrackLookupWriter.hpp"
#include "TrackingPipeline/Io/RootIPPositronTrackParametersReader.hpp"
#include "TrackingPipeline/Io/RootTrackLookupValidationWriter.hpp"
#include "TrackingPipeline/MagneticField/CompositeMagField.hpp"
#include "TrackingPipeline/MagneticField/DipoleMagField.hpp"
#include "TrackingPipeline/MagneticField/QuadrupoleMagField.hpp"
#include "TrackingPipeline/Simulation/E320IPPositronGenerator.hpp"
#include "TrackingPipeline/Simulation/IdealDigitizer.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"
#include "TrackingPipeline/Simulation/RangedUniformMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/StationaryVertexGenerator.hpp"
#include "TrackingPipeline/TrackFinding/TrackLookupEstimationAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/TrackLookupProvider.hpp"
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
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_sim/"
      "E320Pipeline_gdmls/"
      "ettgeom_magnet_pdc_tracker.gdml";
  std::vector<std::string> names{"OPPPSensitive", "DetChamberWindow"};

  Acts::GeometryIdentifier pdcWindowId;
  pdcWindowId.setApproach(1);
  std::vector<Acts::GeometryIdentifier> materialVeto{pdcWindowId};

  std::string materialPath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_sim/"
      "E320Pipeline_material/"
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
  // Cluster filter setup
  SimpleSourceLink::SurfaceAccessor surfaceAccessor{detector.get()};
  auto hourglassFilter = std::make_shared<HourglassFilter>();
  hourglassFilter->surfaceAccessor = surfaceAccessor;

  // --------------------------------------------------------------
  // Event reading

  // Setup the sequencer
  Sequencer::Config seqCfg;
  // seqCfg.events = 2e5;
  seqCfg.numThreads = 1;
  seqCfg.trackFpes = false;
  Sequencer sequencer(seqCfg);

  // Add the sim data reader
  E320Io::E320RootSimDataReader::Config readerCfg = E320Io::defaultSimConfig();
  readerCfg.clusterFilter = hourglassFilter;
  readerCfg.outputSourceLinks = "SimMeasurements";
  readerCfg.outputSimClusters = "SimClusters";
  std::string pathToDir =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_sim/"
      "E320Pipeline_dataInRootFormat/temp";

  // Get the paths to the files in the directory
  for (const auto& entry :
       std::filesystem::recursive_directory_iterator(pathToDir)) {
    if (entry.path().extension() != ".root") {
      continue;
    }
    std::string pathToFile = entry.path();
    readerCfg.filePaths.push_back(pathToFile);
  }

  // The events are not sorted in the directory
  // but we need to process them in order
  std::ranges::sort(readerCfg.filePaths, [](const std::string& a,
                                            const std::string& b) {
    std::size_t idxRootA = a.find_last_of('.');
    std::size_t idxEventA = a.find_last_of('t', idxRootA);
    std::string eventSubstrA = a.substr(idxEventA + 1, idxRootA - idxEventA);

    std::size_t idxRootB = b.find_last_of('.');
    std::size_t idxEventB = b.find_last_of('t', idxRootB);
    std::string eventSubstrB = b.substr(idxEventB + 1, idxRootB - idxEventB);

    return std::stoul(eventSubstrA) < std::stoul(eventSubstrB);
  });

  // Add the reader to the sequencer
  sequencer.addReader(
      std::make_shared<E320Io::E320RootSimDataReader>(readerCfg, logLevel));

  /*// Setup dummy reader*/
  /*DummyReader::Config readerCfg;*/
  /*readerCfg.outputSimClusters = "SimClusters";*/
  /*readerCfg.outputSourceLinks = "SimMeasurements";*/
  /*readerCfg.nEvents = 1e4;*/
  /**/
  /*sequencer.addReader(std::make_shared<DummyReader>(readerCfg));*/

  //  // Setup the measurements creator
  //  Acts::Experimental::DetectorNavigator::Config navCfg;
  //  navCfg.detector = detector.get();
  //  navCfg.resolvePassive = false;
  //  navCfg.resolveMaterial = true;
  //  navCfg.resolveSensitive = true;
  //
  //  Acts::Experimental::DetectorNavigator navigator(
  //      navCfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));
  //  Acts::EigenStepper<> stepper(field);
  //
  //  auto propagator = Propagator(std::move(stepper), std::move(navigator));
  //
  //  auto momGen = std::make_shared<RangedUniformMomentumGenerator>();
  //  momGen->Pranges = {{0.5_GeV, 1.0_GeV}, {1.0_GeV, 1.5_GeV},
  //  {1.5_GeV, 2.0_GeV},
  //                     {2.0_GeV, 2.5_GeV}, {2.5_GeV, 3.0_GeV},
  //                     {3.0_GeV, 3.5_GeV}, {3.5_GeV, 4.0_GeV},
  //                     {4.0_GeV, 4.5_GeV}};
  //
  //  auto digiizer = std::make_shared<IdealDigitizer>();
  //
  //  MeasurementsCreator::Config mcCfg;
  //  mcCfg.vertexGenerator = std::make_shared<StationaryVertexGenerator>();
  //  mcCfg.momentumGenerator = momGen;
  //  mcCfg.hitDigitizer = digiizer;
  //  mcCfg.maxSteps = 300;
  //  mcCfg.isSignal = true;
  //
  //  auto measurementsCreator =
  //      std::make_shared<MeasurementsCreator>(propagator, mcCfg);
  //
  //  MeasurementsEmbeddingAlgorithm::Config mcaCfg;
  //  mcaCfg.inputSourceLinks = "SimMeasurements";
  //  mcaCfg.inputSimClusters = "SimClusters";
  //  mcaCfg.outputSourceLinks = "Measurements";
  //  mcaCfg.outputSimClusters = "Clusters";
  //  mcaCfg.measurementGenerator = measurementsCreator;
  //  mcaCfg.randomNumberSvc =
  //      std::make_shared<RandomNumbers>(RandomNumbers::Config());
  //  mcaCfg.nMeasurements = 1;
  //
  //  sequencer.addAlgorithm(
  //      std::make_shared<MeasurementsEmbeddingAlgorithm>(mcaCfg, logLevel));

  // --------------------------------------------------------------
  // Lookup data validation

  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> refLayers;
  const auto& refVolume = detector->findDetectorVolume("layer0");
  for (const auto* surf : refVolume->surfaces()) {
    refLayers.try_emplace(surf->geometryId(), surf);
  }

  JsonTrackLookupReader::Config lookupReaderCfg;
  lookupReaderCfg.refLayers = refLayers;
  lookupReaderCfg.bins = {1000, 5};

  // Validation algorithm
  TrackLookupValidationAlgorithm::Config validatorCfg;

  TrackLookupProvider::Config providerCfg;
  providerCfg.lookupPath = "lookup-parmigan-1000x5.json";
  providerCfg.trackLookupReader =
      std::make_shared<JsonTrackLookupReader>(lookupReaderCfg);
  TrackLookupProvider provider(providerCfg);

  validatorCfg.refLayers = refLayers;
  validatorCfg.estimator.connect<&TrackLookupProvider::lookup>(&provider);
  validatorCfg.inputClusters = "SimClusters";
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
  validationWriterCfg.path = "lookup-validation-data-1000x5.root";
  validationWriterCfg.treeName = "validation";

  sequencer.addWriter(std::make_shared<RootTrackLookupValidationWriter>(
      validationWriterCfg, logLevel));

  return sequencer.run();
}
