#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Plugins/Json/JsonMaterialDecorator.hpp"
#include "Acts/Plugins/Json/MaterialMapJsonConverter.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/PdgParticle.hpp>
#include <Acts/Material/BinnedSurfaceMaterial.hpp>
#include <Acts/Material/HomogeneousSurfaceMaterial.hpp>
#include <Acts/Material/Interactions.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/RectangleBounds.hpp>

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Alignment/detail/AlignmentStoreBuilders.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryOptions.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Geometry/detail/SurfaceParameters.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Infrastructure/TypeDefinitions.hpp"
#include "TrackingPipeline/Io/DummyReader.hpp"
#include "TrackingPipeline/Io/E320MagneticFieldWriter.hpp"
#include "TrackingPipeline/Io/E320RootSimClusterWriter.hpp"
#include "TrackingPipeline/Simulation/ClusterSizeBasedDigitizer.hpp"
#include "TrackingPipeline/Simulation/ConstantMultiplicityGenerator.hpp"
#include "TrackingPipeline/Simulation/ConvergingBeamGenerator.hpp"
#include "TrackingPipeline/Simulation/ExtendedSourceLinkCreator.hpp"
#include "TrackingPipeline/Simulation/GaussianVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"
#include "TrackingPipeline/Simulation/SimpleDigitizer.hpp"
#include "TrackingPipeline/Simulation/SimpleSourceLinkCreator.hpp"
#include "TrackingPipeline/Simulation/SphericalMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/SurfaceRangedDigitizer.hpp"
#include "TrackingPipeline/Simulation/UniformBackgroundCreator.hpp"
#include "toml++/toml.hpp"

using namespace Acts::UnitLiterals;

std::unique_ptr<const E320::GeometryOptions> E320::GeometryOptions::m_instance =
    nullptr;

int main() {
  const auto& goInst = *E320::GeometryOptions::instance();

  const std::string pathToCfg =
      "/home/romanurmanov/work/TrackingPipeline/ActsTrackingPipeline/conf/"
      "runs/"
      "FastSimRunConfig.toml";
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

  // Contexts
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
      std::cout << surf->polyhedronRepresentation(gctx, 1000).extent() << "\n";
    }
  }

  std::vector<const Acts::Surface*> detSurfaces;
  for (const auto* vol : detector->volumes()) {
    for (const auto* surf : vol->surfaces()) {
      if (surf->geometryId().sensitive() != 0u) {
        detSurfaces.push_back(surf);
      }
    }
  }
  std::vector<SurfaceParameters> constraintSurfaceParameters{
      goInst.ipSurfaceParameters, goInst.beWindowParameters,
      goInst.bpm0Parameters,      goInst.bpm1Parameters,
      goInst.bpm2Parameters,      goInst.bpm3Parameters};

  // --------------------------------------------------------------
  // Alignment setup
  Acts::Vector3 globalShifts(0,
                             getEntryDouble("Geometry", "globalShiftY") * 1_mm,
                             getEntryDouble("Geometry", "globalShiftZ") * 1_mm);
  std::unordered_map<int, Acts::Vector3> localShifts{
      {10,
       Acts::Vector3(0_mm, getEntryDouble("Geometry", "localShiftY10") * 1_um,
                     getEntryDouble("Geometry", "localShiftZ10") * 1_um)},
      {12,
       Acts::Vector3(0_mm, getEntryDouble("Geometry", "localShiftY12") * 1_um,
                     getEntryDouble("Geometry", "localShiftZ12") * 1_um)},
      {14,
       Acts::Vector3(0_mm, getEntryDouble("Geometry", "localShiftY14") * 1_um,
                     getEntryDouble("Geometry", "localShiftZ14") * 1_um)},
      {16,
       Acts::Vector3(0_mm, getEntryDouble("Geometry", "localShiftY16") * 1_um,
                     getEntryDouble("Geometry", "localShiftZ16") * 1_um)},
      {18,
       Acts::Vector3(0_mm, getEntryDouble("Geometry", "localShiftY18") * 1_um,
                     getEntryDouble("Geometry", "localShiftZ18") * 1_um)}};
  Acts::Vector3 globalAngles(
      0_rad, 0_rad, getEntryDouble("Geometry", "globalAngleZ") * 1_mrad);
  std::unordered_map<int, Acts::Vector3> localAngles{
      {10, Acts::Vector3(0_rad, 0_rad,
                         getEntryDouble("Geometry", "localAngleZ10") * 1_mrad)},
      {12, Acts::Vector3(0_rad, 0_rad,
                         getEntryDouble("Geometry", "localAngleZ12") * 1_mrad)},
      {14, Acts::Vector3(0_rad, 0_rad,
                         getEntryDouble("Geometry", "localAngleZ14") * 1_mrad)},
      {16, Acts::Vector3(0_rad, 0_rad,
                         getEntryDouble("Geometry", "localAngleZ16") * 1_mrad)},
      {18,
       Acts::Vector3(0_rad, 0_rad,
                     getEntryDouble("Geometry", "localAngleZ18") * 1_mrad)}};

  auto aStore =
      detail::makeAlignmentStore(gctx, detector.get(), globalShifts,
                                 localShifts, globalAngles, localAngles);
  AlignmentContext alignCtx(aStore);

  // Print alignment parameters
  Acts::GeometryContext defaultGctx;
  gctx = Acts::GeometryContext{alignCtx};
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive() != 0u) {
        std::cout << "-----------------------------------\n";
        std::cout << "SURFACE " << s->geometryId() << "\n";
        std::cout << "CENTER " << s->center(gctx).transpose() << " -- "
                  << s->center(defaultGctx).transpose() << "\n";
        std::cout << "DELTA "
                  << (s->center(gctx) - s->center(defaultGctx)).transpose() *
                         1e3
                  << "\n";
        std::cout << "NORMAL "
                  << s->normal(gctx, s->center(gctx), Acts::Vector3::UnitY())
                         .transpose()
                  << " -- "
                  << s->normal(gctx, s->center(defaultGctx),
                               Acts::Vector3::UnitY())
                         .transpose()
                  << "\n";
        std::cout << "ROTATION \n"
                  << s->transform(gctx).rotation() << " -- \n"
                  << "\n"
                  << s->transform(defaultGctx).rotation() << "\n";

        std::cout << "EXTENT "
                  << s->polyhedronRepresentation(gctx, 1000).extent()
                  << "\n -- \n"
                  << s->polyhedronRepresentation(defaultGctx, 1000).extent()
                  << "\n";
      }
    }
  }

  // --------------------------------------------------------------
  // The magnetic field setup

  auto field = E320::buildMagField(gctx);

  // --------------------------------------------------------------
  // Reference surface
  double halfX = std::numeric_limits<double>::max();
  double halfY = std::numeric_limits<double>::max();

  Acts::RotationMatrix3 refSurfToWorldRotationX =
      Acts::AngleAxis3(goInst.toWorldAngleX, Acts::Vector3::UnitX())
          .toRotationMatrix();
  Acts::RotationMatrix3 refSurfToWorldRotationY =
      Acts::AngleAxis3(goInst.toWorldAngleY, Acts::Vector3::UnitY())
          .toRotationMatrix();
  Acts::RotationMatrix3 refSurfToWorldRotationZ =
      Acts::AngleAxis3(goInst.toWorldAngleZ, Acts::Vector3::UnitZ())
          .toRotationMatrix();

  Acts::Transform3 refSurfTransform = Acts::Transform3::Identity();
  refSurfTransform.translation() = Acts::Vector3::Zero();
  refSurfTransform.rotate(refSurfToWorldRotationX);
  refSurfTransform.rotate(refSurfToWorldRotationY);
  refSurfTransform.rotate(refSurfToWorldRotationZ);

  auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      refSurfTransform, std::make_shared<Acts::RectangleBounds>(halfX, halfY));

  Acts::GeometryIdentifier geoId;
  geoId.setExtra(1);
  refSurface->assignGeometryId(geoId);

  // --------------------------------------------------------------
  // Event reading
  SimpleSourceLink::SurfaceAccessor surfaceAccessor{detector.get()};

  // Setup the sequencer
  Sequencer::Config seqCfg;
  seqCfg.skip = getEntrySizeT("Sequencer", "skip");
  seqCfg.events = getEntrySizeT("Sequencer", "events");
  seqCfg.numThreads = getEntrySizeT("Sequencer", "numThreads");
  seqCfg.trackFpes = getEntryBool("Sequencer", "trackFpes");
  seqCfg.logLevel = logLevel;
  Sequencer sequencer(seqCfg);

  sequencer.addContextDecorator(
      std::make_shared<GeometryContextDecorator>(aStore));

  // --------------------------------------------------------------
  // Add dummy reader
  DummyReader::Config dummyReaderCfg;
  dummyReaderCfg.outputSourceLinks =
      getEntryStr("DummyReader", "outputSourceLinks");
  dummyReaderCfg.outputSimClusters =
      getEntryStr("DummyReader", "outputSimClusters");
  dummyReaderCfg.outputSourceLinkIndices =
      getEntryStr("DummyReader", "outputSourceLinkIndices");
  dummyReaderCfg.nEvents = getEntrySizeT("DummyReader", "nEvents");

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

  // Surface-specific track hit digitizer
  // SurfaceRangedDigitizer::Config trackHitDigitizerCfg;
  // for (const auto& pars : constraintSurfaceParameters) {
  //   trackHitDigitizerCfg.resolutions.insert(
  //       {pars.geoId,
  //        {getEntryDouble("TrackHitDigitizer", "constraintSurfaceResolutionX")
  //        *
  //             1_um,
  //         getEntryDouble("TrackHitDigitizer", "constraintSurfaceResolutionY")
  //         *
  //             1_um}});
  // }
  // for (const auto& pars : goInst.tcParameters) {
  //   trackHitDigitizerCfg.resolutions.insert(
  //       {pars.geoId,
  //        {getEntryDouble("TrackHitDigitizer", "trackingSurfaceResolutionX") *
  //             1_um,
  //         getEntryDouble("TrackHitDigitizer", "trackingSurfaceResolutionY") *
  //             1_um}});
  // }
  // auto trackHitDigitizer =
  //     std::make_shared<SurfaceRangedDigitizer>(trackHitDigitizerCfg);

  // Track angle digitizer
  SimpleDigitizer::Config angleDigitizerCfg;
  angleDigitizerCfg.resolution = {
      getEntryDouble("TrackAngleDigitizer", "constraintSurfaceResolutionPhi") *
          1_mrad,
      getEntryDouble("TrackAngleDigitizer",
                     "constraintSurfaceResolutionTheta") *
          1_mrad};
  auto angleDigitizer = std::make_shared<SimpleDigitizer>(angleDigitizerCfg);

  // Cluster-size-dependent track hit digitizer
  ClusterSizeBasedDigitizer::Config trackHitDigitizerCfg;
  const auto* clusterSizeDigitizationPars =
      runCfg["ClusterSizeBasedDigitizer"].as_array();

  for (auto it = clusterSizeDigitizationPars->begin();
       it != clusterSizeDigitizationPars->end(); it++) {
    const auto& entry = *it->as_table();
    trackHitDigitizerCfg.clSizeProbsStdDevs.insert(
        {entry["clusterSize"].value<std::size_t>().value(),
         {
             entry["samplingProb"].value<double>().value(),
             entry["resolutionX"].value<double>().value() * 1_um,
             entry["resolutionY"].value<double>().value() * 1_um,
         }});
  }
  auto trackHitDigitizer =
      std::make_shared<ClusterSizeBasedDigitizer>(trackHitDigitizerCfg);

  // Vertex generator
  double vertexRes =
      getEntryDouble("GaussianVertexGenerator", "vertexResolution") * 1_um;

  GaussianVertexGenerator::Config vertexGenCfg;
  // vertexGenCfg.mean = Acts::Vector3(goInst.bpm0CenterPrimary - 5_mm, 0, 0);
  vertexGenCfg.mean = Acts::Vector3(0, 0, 0);
  // Acts::Vector3(goInst.ipSurfaceCenterPrimary - 0.1_mm, goInst.tcCenterLong,
  //               goInst.tcCenterShort);
  vertexGenCfg.cov = Acts::SquareMatrix3::Identity() * vertexRes * vertexRes;
  auto vertexGen = std::make_shared<GaussianVertexGenerator>(vertexGenCfg);

  // Momentum generator
  double minAbsP =
      getEntryDouble("SphericalMomentumGenerator", "minAbsoluteMomentum") *
      1_GeV;
  double maxAbsP =
      getEntryDouble("SphericalMomentumGenerator", "maxAbsoluteMomentum") *
      1_GeV;

  double minPhi =
      getEntryDouble("SphericalMomentumGenerator", "minPhi") * 1_mrad;
  double maxPhi =
      getEntryDouble("SphericalMomentumGenerator", "maxPhi") * 1_mrad;

  double minTheta =
      getEntryDouble("SphericalMomentumGenerator", "minTheta") * 1_mrad;
  double maxTheta =
      getEntryDouble("SphericalMomentumGenerator", "maxTheta") * 1_mrad;

  SphericalMomentumGenerator::Config momGenCfg;
  momGenCfg.pRange = {minAbsP, maxAbsP};
  momGenCfg.phiRange = {minPhi, maxPhi};
  momGenCfg.thetaRange = {M_PI_2 + minTheta, M_PI_2 + maxTheta};

  auto momGen = std::make_shared<SphericalMomentumGenerator>(momGenCfg);

  // // Beam generator
  // ConvergingBeamGenerator::Config beamGenCfg;
  // beamGenCfg.primaryIdx = goInst.primaryIdx;
  // beamGenCfg.longIdx = goInst.longIdx;
  // beamGenCfg.shortIdx = goInst.shortIdx;
  // beamGenCfg.referencePositionPrimary = goInst.bpm0CenterPrimary - 5_mm;
  // beamGenCfg.waistPosition =
  //     Acts::Vector3(goInst.ipSurfaceCenterPrimary - 0.1_mm, 0, 0_mm);
  // beamGenCfg.waistSigmaLong = 30_um;
  // beamGenCfg.waistSigmaShort = 30_um;
  // beamGenCfg.waistMeanThetaLong = 0_rad;
  // beamGenCfg.waistMeanThetaShort = 0_rad;
  // beamGenCfg.waistSigmaThetaLong = 1e-3_rad;
  // beamGenCfg.waistSigmaThetaShort = 1e-3_rad;
  // // beamGenCfg.beamEnergyMin = 3.0_GeV;
  // // beamGenCfg.beamEnergyMax = 4.0_GeV;
  // beamGenCfg.beamEnergyMin = 1.5_GeV;
  // beamGenCfg.beamEnergyMax = 2.5_GeV;
  // auto beamGen = std::make_shared<ConvergingBeamGenerator>(beamGenCfg);

  // // Beam positrons generator
  // RootVertexMomentumReaderGenerator::Config beamGenCfg;
  // beamGenCfg.filePaths = {
  //     "/home/romanurmanov/work/E320/E320Prototype/"
  //     "E320Prototype_dataInRootFormat/sim/alignment/global/"
  //     "ip_data_aligned_beamline/fitted-beam-sim-test.root"};
  // beamGenCfg.treeName = "fitted-beam";
  // beamGenCfg.vertexBranch = "originPositionTruth";
  // beamGenCfg.directionBranch = "originMomentumTruth";
  // beamGenCfg.absMomentumMin = 3.0_GeV;
  // beamGenCfg.absMomentumMax = 4.0_GeV;
  // beamGenCfg.startIdx = 5 * 2.5e4;
  // auto beamGen =
  //     std::make_shared<RootVertexMomentumReaderGenerator>(beamGenCfg);

  // Source link creators
  SimpleSourceLinkCreator::Config simpleSourceLinkCreatorCfg;
  simpleSourceLinkCreatorCfg.hitDigitizer = trackHitDigitizer;
  SimpleSourceLinkCreator simpleSourceLinkCreator(simpleSourceLinkCreatorCfg);

  ExtendedSourceLinkCreator::Config extendedSourceLinkCreatorCfg;
  extendedSourceLinkCreatorCfg.hitDigitizer = trackHitDigitizer;
  extendedSourceLinkCreatorCfg.angleDigitizer = angleDigitizer;
  ExtendedSourceLinkCreator extendedSourceLinkCreator(
      extendedSourceLinkCreatorCfg);

  // Measurement creator
  MeasurementsCreator::Config measCreatorCfg{
      .vertexGenerator = vertexGen,
      .momentumGenerator = momGen,
      .referenceSurface = refSurface.get(),
      .maxSteps = 1000,
      .isSignal = true,
      .hypothesis = Acts::ParticleHypothesis::electron(),
      .charge = 1_e};

  for (const auto& pars : constraintSurfaceParameters) {
    if (pars.geoId != goInst.ipSurfaceParameters.geoId) {
      measCreatorCfg.sourceLinkCreators[pars.geoId]
          .connect<&SimpleSourceLinkCreator::operator()>(
              &simpleSourceLinkCreator);
    } else {
      measCreatorCfg.sourceLinkCreators[pars.geoId]
          .connect<&ExtendedSourceLinkCreator::operator()>(
              &extendedSourceLinkCreator);
    }
  }
  for (const auto& pars : goInst.tcParameters) {
    measCreatorCfg.sourceLinkCreators[pars.geoId]
        .connect<&SimpleSourceLinkCreator::operator()>(
            &simpleSourceLinkCreator);
  }

  std::unordered_map<Acts::GeometryIdentifier, MeasurementsCreator::Constraints>
      measCreatorConstraints;
  for (auto& det : detector->detectorElements()) {
    const auto& surface = det->surface();
    const auto& geoId = surface.geometryId();
    if (geoId.sensitive() >= goInst.bpm0Parameters.geoId &&
        geoId.sensitive() <= goInst.bpm3Parameters.geoId) {
      measCreatorConstraints.insert({geoId, {-3e1, 3e1, -3e1, 3e1}});
    }
  }
  measCreatorCfg.constraints = measCreatorConstraints;

  auto measCreator = std::make_shared<MeasurementsCreator>(
      measCreatorPropagator, measCreatorCfg);

  /// Event multiplicity generator
  ConstantMultiplicityGenerator::Config sigMultiplicityGeneratorCfg{};
  sigMultiplicityGeneratorCfg.eventMultiplicity =
      getEntrySizeT("SigConstantMultiplicityGenerator", "nMeasurements");

  auto sigMultiplicityGenerator =
      std::make_shared<ConstantMultiplicityGenerator>(
          sigMultiplicityGeneratorCfg);

  /// Measurement embedding algorithm
  MeasurementsEmbeddingAlgorithm::Config measCreatorAlgoCfg;
  measCreatorAlgoCfg.inputSourceLinks =
      getEntryStr("SigMeasurementsEmbeddingAlgorithm", "inputSourceLinks");
  measCreatorAlgoCfg.inputSimClusters =
      getEntryStr("SigMeasurementsEmbeddingAlgorithm", "inputSimClusters");
  measCreatorAlgoCfg.inputSourceLinkIndices = getEntryStr(
      "SigMeasurementsEmbeddingAlgorithm", "inputSourceLinkIndices");
  measCreatorAlgoCfg.outputSourceLinks =
      getEntryStr("SigMeasurementsEmbeddingAlgorithm", "outputSourceLinks");
  measCreatorAlgoCfg.outputSimClusters =
      getEntryStr("SigMeasurementsEmbeddingAlgorithm", "outputSimClusters");
  measCreatorAlgoCfg.outputSourceLinkIndices = getEntryStr(
      "SigMeasurementsEmbeddingAlgorithm", "outputSourceLinkIndices");
  measCreatorAlgoCfg.measurementGenerator = measCreator;
  measCreatorAlgoCfg.multiplicityGenerator = sigMultiplicityGenerator;
  measCreatorAlgoCfg.randomNumberSvc =
      std::make_shared<RandomNumbers>(RandomNumbers::Config());

  sequencer.addAlgorithm(std::make_shared<MeasurementsEmbeddingAlgorithm>(
      measCreatorAlgoCfg, logLevel));

  // --------------------------------------------------------------
  // Add background

  // Background creator
  UniformBackgroundCreator::Config bkgCreatorCfg;
  bkgCreatorCfg.resolution = {
      getEntryDouble("UniformBackgroundCreator", "resolutionX") * 1_um,
      getEntryDouble("UniformBackgroundCreator", "resolutionY") * 1_um};
  bkgCreatorCfg.nMeasurementsPerSurface =
      getEntrySizeT("UniformBackgroundCreator", "nMeasurementsPerSurface");
  bkgCreatorCfg.surfaces = detSurfaces;

  auto bkgCreator = std::make_shared<UniformBackgroundCreator>(bkgCreatorCfg);

  /// Event multiplicity generator
  ConstantMultiplicityGenerator::Config bkgMultiplicityGeneratorCfg{};
  bkgMultiplicityGeneratorCfg.eventMultiplicity =
      getEntrySizeT("BkgConstantMultiplicityGenerator", "nMeasurements");

  auto bkgMultiplicityGenerator =
      std::make_shared<ConstantMultiplicityGenerator>(
          bkgMultiplicityGeneratorCfg);

  MeasurementsEmbeddingAlgorithm::Config bkgCreatorAlgoCfg;
  bkgCreatorAlgoCfg.inputSourceLinks =
      getEntryStr("BkgMeasurementsEmbeddingAlgorithm", "inputSourceLinks");
  bkgCreatorAlgoCfg.inputSimClusters =
      getEntryStr("BkgMeasurementsEmbeddingAlgorithm", "inputSimClusters");
  bkgCreatorAlgoCfg.inputSourceLinkIndices = getEntryStr(
      "BkgMeasurementsEmbeddingAlgorithm", "inputSourceLinkIndices");
  bkgCreatorAlgoCfg.outputSourceLinks =
      getEntryStr("BkgMeasurementsEmbeddingAlgorithm", "outputSourceLinks");
  bkgCreatorAlgoCfg.outputSimClusters =
      getEntryStr("BkgMeasurementsEmbeddingAlgorithm", "outputSimClusters");
  bkgCreatorAlgoCfg.outputSourceLinkIndices = getEntryStr(
      "BkgMeasurementsEmbeddingAlgorithm", "outputSourceLinkIndices");
  bkgCreatorAlgoCfg.measurementGenerator = bkgCreator;
  measCreatorAlgoCfg.multiplicityGenerator = bkgMultiplicityGenerator;
  bkgCreatorAlgoCfg.randomNumberSvc =
      std::make_shared<RandomNumbers>(RandomNumbers::Config());

  // sequencer.addAlgorithm(std::make_shared<MeasurementsEmbeddingAlgorithm>(
  //     bkgCreatorAlgoCfg, logLevel));

  // --------------------------------------------------------------
  // Event write out

  // Sim cluster writer
  E320::E320RootSimClusterWriter::Config clusterWriterCfgSig;
  clusterWriterCfgSig.inputClusters =
      getEntryStr("E320RootSimClusterWriter", "inputClusters");
  clusterWriterCfgSig.treeName =
      getEntryStr("E320RootSimClusterWriter", "treeName");
  clusterWriterCfgSig.filePath =
      getEntryStr("E320RootSimClusterWriter", "filePath");

  sequencer.addWriter(std::make_shared<E320::E320RootSimClusterWriter>(
      clusterWriterCfgSig, logLevel));

  // Magnetic field writer
  E320::E320MagneticFieldWriter::Config magFieldWriterCfg;
  magFieldWriterCfg.treeName =
      getEntryStr("E320MagneticFieldWriter", "treeName");
  magFieldWriterCfg.filePath =
      getEntryStr("E320MagneticFieldWriter", "filePath");

  sequencer.addWriter(std::make_shared<E320::E320MagneticFieldWriter>(
      magFieldWriterCfg, logLevel));

  return sequencer.run();
}
