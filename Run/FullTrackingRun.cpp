#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Seeding/PathSeeder.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <filesystem>
#include <unordered_map>

#include "TrackingPipeline/Clustering/HourglassFilter.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/DummyReader.hpp"
#include "TrackingPipeline/Io/E320RootDataReader.hpp"
#include "TrackingPipeline/Io/JsonTrackLookupReader.hpp"
#include "TrackingPipeline/Io/RootFittedSimTrackWriter.hpp"
#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"
#include "TrackingPipeline/Io/RootSimSeedWriter.hpp"
#include "TrackingPipeline/Io/RootSimTrackCandidateWriter.hpp"
#include "TrackingPipeline/Io/RootTrackParamsReader.hpp"
#include "TrackingPipeline/MagneticField/CompositeMagField.hpp"
#include "TrackingPipeline/MagneticField/DipoleMagField.hpp"
#include "TrackingPipeline/MagneticField/QuadrupoleMagField.hpp"
#include "TrackingPipeline/Simulation/E320CptBkgGenerator.hpp"
#include "TrackingPipeline/Simulation/E320HistDigitizer.hpp"
#include "TrackingPipeline/Simulation/KDEMomentumxGenerator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"
#include "TrackingPipeline/Simulation/PowerLawVertexGenerator.hpp"
#include "TrackingPipeline/TrackFinding/CKFTrackFindingAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/ForwardOrderedIntersectionFinder.hpp"
#include "TrackingPipeline/TrackFinding/LayerPathWidthProvider.hpp"
#include "TrackingPipeline/TrackFinding/PathSeedingAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/SourceLinkGridConstructor.hpp"
#include "TrackingPipeline/TrackFinding/TrackLookupProvider.hpp"
#include "TrackingPipeline/TrackFitting/TrackFittingAlgorithm.hpp"

using ActionList = Acts::ActionList<>;
using AbortList = Acts::AbortList<Acts::EndOfWorldReached>;

using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                    Acts::Experimental::DetectorNavigator>;
using PropagatorOptions =
    typename Propagator::template Options<ActionList, AbortList>;

using Trajectory = Acts::VectorMultiTrajectory;
using KFTrackContainer = Acts::VectorTrackContainer;
using KF = Acts::KalmanFitter<Propagator, Trajectory>;

using CKFTrackContainer = Acts::TrackContainer<Acts::VectorTrackContainer,
                                               Acts::VectorMultiTrajectory,
                                               Acts::detail::ValueHolder>;

using TrackStateContainerBackend =
    typename CKFTrackContainer::TrackStateContainerBackend;

using namespace Acts::UnitLiterals;

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

  Acts::ActsScalar dipoleB = 0.31_T;
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
  hourglassFilter->surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  // --------------------------------------------------------------
  // Event reading

  // Setup the sequencer
  Sequencer::Config seqCfg;
  //   seqCfg.events = 100;
  seqCfg.numThreads = 1;
  seqCfg.trackFpes = false;
  Sequencer sequencer(seqCfg);

  //   // Add the sim data reader
  //   E320Io::E320RootSimDataReader::Config readerCfg =
  //   E320Io::defaultSimConfig(); readerCfg.clusterFilter = hourglassFilter;
  //   readerCfg.outputSourceLinks = "Measurements";
  //   readerCfg.outputSimClusters = "Clusters";
  //   std::string pathToDir =
  //   "/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_dataInRootFormat/"
  //   "temp";

  //   // Get the paths to the files in the directory
  //   for (const auto& entry : std::filesystem::directory_iterator(pathToDir))
  //   {
  // std::string pathToFile = entry.path();
  // readerCfg.filePaths.push_back(pathToFile);
  //   }

  //   // The events are not sorted in the directory
  //   // but we need to process them in order
  //   std::ranges::sort(readerCfg.filePaths, [](const std::string& a,
  // const std::string& b) {
  // std::size_t idxRootA = a.find_last_of('.');
  // std::size_t idxEventA = a.find_last_of('t', idxRootA);
  // std::string eventSubstrA = a.substr(idxEventA + 1, idxRootA - idxEventA);

  // std::size_t idxRootB = b.find_last_of('.');
  // std::size_t idxEventB = b.find_last_of('t', idxRootB);
  // std::string eventSubstrB = b.substr(idxEventB + 1, idxRootB - idxEventB);

  // return std::stoul(eventSubstrA) < std::stoul(eventSubstrB);
  //   });

  //   // Add the reader to the sequencer
  //   sequencer.addReader(
  //   std::make_shared<E320Io::E320RootSimDataReader>(readerCfg, logLevel));

  // --------------------------------------------------------------
  // Add dummy reader
  DummyReader::Config dummyReaderCfg;
  dummyReaderCfg.outputSourceLinks = "SimMeasurements";
  dummyReaderCfg.outputSimClusters = "SimClusters";
  dummyReaderCfg.nEvents = 10;

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

  Acts::ActsScalar y0 = gOpt.chipY.at(0) - gOpt.chipSizeY / 2;
  Acts::ActsScalar y1 = gOpt.chipY.at(8) + gOpt.chipSizeY / 2;

  Acts::ActsScalar x0 = gOpt.chipX - gOpt.chipSizeX / 2;
  Acts::ActsScalar x1 = gOpt.chipX + gOpt.chipSizeX / 2;

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
  cptBkgGenCfg.yPower = -77.0091;
  cptBkgGenCfg.yShift = -10630.9;
  cptBkgGenCfg.zProbs = {0.29, 0.26, 0.24, 0.21};
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
  cptMceCfg.outputSourceLinks = "MeasurementsWithNCS";
  cptMceCfg.outputSimClusters = "ClustersWithNCS";
  cptMceCfg.measurementGenerator = cptMeasurementsCreator;
  cptMceCfg.clusterFilter = hourglassFilter;
  cptMceCfg.randomNumberSvc =
      std::make_shared<RandomNumbers>(RandomNumbers::Config());
  cptMceCfg.nMeasurements = 6600;

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
  auto beamVertexGen = std::make_shared<E320Sim::E320PowerLawVertexGenerator>();

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
  beamMceCfg.inputSourceLinks = "MeasurementsWithNCS";
  beamMceCfg.inputSimClusters = "ClustersWithNCS";
  beamMceCfg.outputSourceLinks = "Measurements";
  beamMceCfg.outputSimClusters = "Clusters";
  beamMceCfg.measurementGenerator = beamMeasurementsCreator;
  beamMceCfg.clusterFilter = hourglassFilter;
  beamMceCfg.randomNumberSvc =
      std::make_shared<RandomNumbers>(RandomNumbers::Config());
  beamMceCfg.nMeasurements = 34500;

  sequencer.addAlgorithm(
      std::make_shared<MeasurementsEmbeddingAlgorithm>(beamMceCfg, logLevel));

  // --------------------------------------------------------------
  // The path seeding setup
  auto pathSeederCfg = Acts::PathSeeder::Config();

  // Combine layers into surfaces to take care of the gaps
  std::vector<std::shared_ptr<Acts::Surface>> layerPtrs;
  for (int i = 0; i < gOpt.staveZ.size(); i++) {
    // Bounds
    double halfX = gOpt.chipSizeX / 2;
    double halfY = ((gOpt.chipY.at(8) + gOpt.chipSizeY / 2) -
                    (gOpt.chipY.at(0) - gOpt.chipSizeY / 2)) /
                   2;

    // Center
    double centerY = ((-gOpt.chipY.at(8) - gOpt.chipSizeY / 2) +
                      (-gOpt.chipY.at(0) + gOpt.chipSizeY / 2)) /
                     2;

    Acts::Transform3 transform(Acts::Translation3(Acts::Vector3(
                                   gOpt.chipX, gOpt.staveZ.at(i), centerY)) *
                               gOpt.actsToWorldRotation.inverse());

    auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
        transform, std::make_shared<Acts::RectangleBounds>(halfX, halfY));

    Acts::GeometryIdentifier geoId;
    geoId.setSensitive(i + 1);
    surface->assignGeometryId(std::move(geoId));
    layerPtrs.push_back(surface);

    if (i == 0) {
      pathSeederCfg.refLayerIds.push_back(geoId);
    }
  }

  // Intersection finder
  std::vector<const Acts::Surface*> layers;
  for (auto layer : layerPtrs) {
    layers.push_back(layer.get());
  }

  ForwardOrderedIntersectionFinder::Config intersectionFinderCfg;
  intersectionFinderCfg.layers = std::move(layers);
  intersectionFinderCfg.tol = (gOpt.chipY.at(1) - gOpt.chipSizeY / 2) -
                              (gOpt.chipY.at(0) + gOpt.chipSizeY / 2) + 1_mm;
  ForwardOrderedIntersectionFinder intersectionFinder(intersectionFinderCfg);

  pathSeederCfg.intersectionFinder
      .connect<&ForwardOrderedIntersectionFinder::operator()>(
          &intersectionFinder);

  // Path width provider
  std::map<std::int32_t, std::pair<Acts::ActsScalar, Acts::ActsScalar>>
      pathWidths = {
          {0, {100_um, 100_um}},
          {1, {300_um, 300_um}},
          {2, {305_um, 350_um}},
          {3, {310_um, 400_um}},
      };

  LayerPathWidthProvider pathWidthProvider(pathWidths);

  pathSeederCfg.pathWidthProvider.connect<&LayerPathWidthProvider::operator()>(
      &pathWidthProvider);

  // Grid to bin the source links
  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> layerMap;
  for (auto layer : layerPtrs) {
    layerMap.insert({layer->geometryId(), layer.get()});
  }

  SourceLinkGridConstructor::Config gridConstructorCfg{
      .bins = std::make_pair(1, 1000), .layers = layerMap};
  gridConstructorCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  auto gridConstructor =
      std::make_shared<SourceLinkGridConstructor>(gridConstructorCfg);

  // Estimator of the IP and first hit
  // parameters of the track
  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> refLayers;
  const auto& refVolume = detector->findDetectorVolume("layer0");
  for (const auto* surf : refVolume->surfaces()) {
    refLayers.try_emplace(surf->geometryId(), surf);
  }

  JsonTrackLookupReader::Config lookupReaderCfg;
  lookupReaderCfg.refLayers = refLayers;
  lookupReaderCfg.bins = {1000, 1};

  TrackLookupProvider::Config lookupProviderCfg;
  lookupProviderCfg.lookupPath =
      "/home/romanurmanov/lab/LUXE/acts_tracking/E320Pipeline_lookups/"
      "RangedUniform_05_45_Stationary_000_1000x1_200k/lookup.json";
  lookupProviderCfg.trackLookupReader =
      std::make_shared<JsonTrackLookupReader>(lookupReaderCfg);
  TrackLookupProvider lookupProvider(lookupProviderCfg);

  pathSeederCfg.trackEstimator.connect<&TrackLookupProvider::lookup>(
      &lookupProvider);

  // Create the path seeder algorithm
  auto seedingAlgoCfg = PathSeedingAlgorithm::Config();
  seedingAlgoCfg.seeder = std::make_shared<Acts::PathSeeder>(pathSeederCfg);
  seedingAlgoCfg.sourceLinkGridConstructor = gridConstructor;
  seedingAlgoCfg.inputSourceLinks = "Measurements";
  seedingAlgoCfg.outputSeeds = "PathSeeds";
  seedingAlgoCfg.minSeedSize = 4;
  seedingAlgoCfg.maxSeedSize = 100;
  sequencer.addAlgorithm(
      std::make_shared<PathSeedingAlgorithm>(seedingAlgoCfg, logLevel));

  // --------------------------------------------------------------
  // Track finding
  Acts::Experimental::DetectorNavigator::Config ckfNavigatorCfg;
  ckfNavigatorCfg.detector = detector.get();
  ckfNavigatorCfg.resolvePassive = false;
  ckfNavigatorCfg.resolveMaterial = true;
  ckfNavigatorCfg.resolveSensitive = true;
  Acts::Experimental::DetectorNavigator ckfNavigator(
      ckfNavigatorCfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));

  Acts::EigenStepper<> ckfStepper(field);
  auto ckfPropagator =
      Propagator(std::move(ckfStepper), std::move(ckfNavigator),
                 Acts::getDefaultLogger("Propagator", logLevel));

  Acts::CombinatorialKalmanFilter<Propagator, CKFTrackContainer> ckf(
      ckfPropagator,
      Acts::getDefaultLogger("CombinatorialKalmanFilter", logLevel));

  // Configuration for the measurement selector
  std::vector<
      std::pair<Acts::GeometryIdentifier, Acts::MeasurementSelectorCuts>>
      cuts;
  for (auto& vol : detector->volumes()) {
    for (auto& surf : vol->surfaces()) {
      if (vol->name() == "layer0") {
        cuts.push_back(
            {surf->geometryId(),
             {{}, {std::numeric_limits<Acts::ActsScalar>::max()}, {1000u}}});
      } else {
        cuts.push_back({surf->geometryId(), {{}, {0.5}, {1u}}});
      }
    }
  }
  Acts::MeasurementSelector::Config measurementSelectorCfg(cuts);

  Acts::MeasurementSelector measSel{measurementSelectorCfg};

  // CKF extensions
  Acts::GainMatrixUpdater ckfUpdater;

  Acts::CombinatorialKalmanFilterExtensions<CKFTrackContainer> ckfExtensions;
  ckfExtensions.calibrator.template connect<
      &simpleSourceLinkCalibrator<TrackStateContainerBackend>>();
  ckfExtensions.updater.template connect<
      &Acts::GainMatrixUpdater::operator()<TrackStateContainerBackend>>(
      &ckfUpdater);
  ckfExtensions.measurementSelector.template connect<
      &Acts::MeasurementSelector::select<TrackStateContainerBackend>>(&measSel);

  CKFTrackFindingAlgorithm<Propagator, CKFTrackContainer>::Config
      trackFindingCfg{
          .ckf = ckf,
      };
  trackFindingCfg.extensions = ckfExtensions;
  trackFindingCfg.inputSeeds = "PathSeeds";
  trackFindingCfg.outputTrackCandidates = "TrackCandidates";
  trackFindingCfg.minCandidateSize = 4;
  trackFindingCfg.maxCandidateSize = 4;

  auto trackFindingAlgorithm =
      std::make_shared<CKFTrackFindingAlgorithm<Propagator, CKFTrackContainer>>(
          trackFindingCfg, logLevel);
  sequencer.addAlgorithm(trackFindingAlgorithm);

  // --------------------------------------------------------------
  // Track fitting

  Acts::GainMatrixUpdater kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;

  // Initialize track fitter options
  Acts::KalmanFitterExtensions<Trajectory> extensions;
  // Add calibrator
  extensions.calibrator.connect<&simpleSourceLinkCalibrator<Trajectory>>();
  // Add the updater
  extensions.updater.connect<&Acts::GainMatrixUpdater::operator()<Trajectory>>(
      &kfUpdater);
  // Add the smoother
  extensions.smoother
      .connect<&Acts::GainMatrixSmoother::operator()<Trajectory>>(&kfSmoother);
  // Add the surface accessor
  extensions.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  auto propOptions = PropagatorOptions(gctx, mctx);

  propOptions.maxSteps = 300;

  auto options =
      Acts::KalmanFitterOptions(gctx, mctx, cctx, extensions, propOptions);

  // Reference surface for sampling the track at the IP
  double halfX = std::numeric_limits<double>::max();
  double halfY = std::numeric_limits<double>::max();

  // double refZ = gOpt.dipoleTranslation.z() + gOpt.dipoleBounds.at(2);
  double refZ = 0;
  Acts::Transform3 transform(Acts::Translation3(Acts::Vector3(0, refZ, 0)) *
                             gOpt.actsToWorldRotation.inverse());

  auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      transform, std::make_shared<Acts::RectangleBounds>(halfX, halfY));

  Acts::GeometryIdentifier geoId;
  geoId.setExtra(1);
  refSurface->assignGeometryId(std::move(geoId));

  options.referenceSurface = refSurface.get();

  Acts::Experimental::DetectorNavigator::Config cfg;
  cfg.detector = detector.get();
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Experimental::DetectorNavigator kfNavigator(
      cfg, Acts::getDefaultLogger("DetectorNavigator", logLevel));

  Acts::EigenStepper<> kfStepper(std::move(field));
  auto kfPropagator =
      Propagator(std::move(kfStepper), std::move(kfNavigator),
                 Acts::getDefaultLogger("Propagator", logLevel));

  const auto fitter = KF(
      kfPropagator, Acts::getDefaultLogger("DetectorKalmanFilter", logLevel));

  // Add the track fitting algorithm to the sequencer
  TrackFittingAlgorithm<Propagator, Trajectory, KFTrackContainer>::Config
      fitterCfg{.inputCollection = "TrackCandidates",
                .outputCollection = "Tracks",
                .fitter = fitter,
                .kfOptions = options};

  sequencer.addAlgorithm(
      std::make_shared<
          TrackFittingAlgorithm<Propagator, Trajectory, KFTrackContainer>>(
          fitterCfg, logLevel));

  // --------------------------------------------------------------
  // Event write out

  // Sim cluster writer
  auto clusterWriterCfg = RootSimClusterWriter::Config();

  clusterWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  clusterWriterCfg.inputClusters = "Clusters";
  clusterWriterCfg.treeName = "clusters";
  clusterWriterCfg.filePath = "clusters-bkg-10-nominal.root";

  sequencer.addWriter(
      std::make_shared<RootSimClusterWriter>(clusterWriterCfg, logLevel));

  // Seed writer
  auto seedWriterCfg = RootSimSeedWriter::Config();

  seedWriterCfg.inputSeeds = "PathSeeds";
  seedWriterCfg.inputTruthClusters = "Clusters";
  seedWriterCfg.treeName = "seeds";
  seedWriterCfg.filePath = "seeds-bkg-10-nominal.root";
  seedWriterCfg.targetTrueTrackSize = 4;

  sequencer.addWriter(
      std::make_shared<RootSimSeedWriter>(seedWriterCfg, logLevel));

  // Track candidate writer
  auto trackCandidateWriterCfg = RootSimTrackCandidateWriter::Config();

  trackCandidateWriterCfg.inputTrackCandidates = "TrackCandidates";
  trackCandidateWriterCfg.inputTruthClusters = "Clusters";
  trackCandidateWriterCfg.treeName = "track-candidates";
  trackCandidateWriterCfg.filePath = "track-candidates-bkg-10-nominal.root";
  trackCandidateWriterCfg.targetTrueTrackSize = 4;

  sequencer.addWriter(std::make_shared<RootSimTrackCandidateWriter>(
      trackCandidateWriterCfg, logLevel));

  auto trackWriterCfg = RootFittedSimTrackWriter::Config();
  trackWriterCfg.surfaceAccessor
      .connect<&SimpleSourceLink::SurfaceAccessor::operator()>(
          &surfaceAccessor);

  trackWriterCfg.inputKFTracks = "Tracks";
  trackWriterCfg.inputTruthClusters = "Clusters";
  trackWriterCfg.treeName = "fitted-tracks";
  trackWriterCfg.filePath = "fitted-tracks-bkg-10-nominal.root";
  trackWriterCfg.targetTrueTrackSize = 4;

  sequencer.addWriter(
      std::make_shared<RootFittedSimTrackWriter>(trackWriterCfg, logLevel));

  return sequencer.run();
}
