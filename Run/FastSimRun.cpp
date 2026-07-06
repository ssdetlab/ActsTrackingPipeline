#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/RectangleBounds.hpp>

#include <memory>
#include <string>
#include <vector>

#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Alignment/detail/AlignmentStoreBuilders.hpp"
#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/GeometryContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Infrastructure/TypeDefinitions.hpp"
#include "TrackingPipeline/Io/DummyReader.hpp"
#include "TrackingPipeline/Io/E320MagneticFieldWriter.hpp"
#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"
#include "TrackingPipeline/Simulation/ClusterSizeBasedDigitizer.hpp"
#include "TrackingPipeline/Simulation/GaussianVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"
#include "TrackingPipeline/Simulation/SimpleDigitizer.hpp"
#include "TrackingPipeline/Simulation/SphericalMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/SurfaceRangedDigitizer.hpp"
#include "TrackingPipeline/Simulation/UniformBackgroundCreator.hpp"

using namespace Acts::UnitLiterals;

std::unique_ptr<const E320::GeometryOptions> E320::GeometryOptions::m_instance =
    nullptr;

int main() {
  const auto& goInst = *E320::GeometryOptions::instance();

  // Set the log level
  Acts::Logging::Level logLevel = Acts::Logging::INFO;

  // Contexts
  Acts::GeometryContext gctx;

  // --------------------------------------------------------------
  // Detector setup

  auto detector = E320::buildDetector(gctx, nullptr);

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

  // --------------------------------------------------------------
  // Alignment setup
  Acts::Vector3 globalShiftMean(0, 0_mm, 0_mm);
  Acts::Vector3 globalShiftStdErr(0, 0_mm, 0_mm);

  std::unordered_map<int, Acts::Vector3> localShiftsMean{
      {10, Acts::Vector3(0_mm, 0_um, 0_um)},
      {12, Acts::Vector3(0_mm, 0_um, 0_um)},
      {14, Acts::Vector3(0_mm, 0_um, 0_um)},
      {16, Acts::Vector3(0_mm, 0_um, 0_um)},
      {18, Acts::Vector3(0_mm, 0_um, 0_um)}};
  std::unordered_map<int, Acts::Vector3> localShiftsStdErr{
      {10, Acts::Vector3(0_mm, 0_um, 0_um)},
      {12, Acts::Vector3(0_mm, 0_um, 0_um)},
      {14, Acts::Vector3(0_mm, 0_um, 0_um)},
      {16, Acts::Vector3(0_mm, 0_um, 0_um)},
      {18, Acts::Vector3(0_mm, 0_um, 0_um)}};

  Acts::Vector3 globalAnglesMean(0_rad, 0_rad, 0_rad);
  Acts::Vector3 globalAnglesStdErr(0_rad, 0_rad, 0_rad);

  std::unordered_map<int, Acts::Vector3> localAnglesMean{
      {10, Acts::Vector3(0_rad, 0_rad, 0_rad)},
      {12, Acts::Vector3(0_rad, 0_rad, 0_rad)},
      {14, Acts::Vector3(0_rad, 0_rad, 0_rad)},
      {16, Acts::Vector3(0_rad, 0_rad, 0_rad)},
      {18, Acts::Vector3(0_rad, 0_rad, 0_rad)}};
  std::unordered_map<int, Acts::Vector3> localAnglesStdErr{
      {10, Acts::Vector3(0_rad, 0_rad, 0_rad)},
      {12, Acts::Vector3(0_rad, 0_rad, 0_rad)},
      {14, Acts::Vector3(0_rad, 0_rad, 0_rad)},
      {16, Acts::Vector3(0_rad, 0_rad, 0_rad)},
      {18, Acts::Vector3(0_rad, 0_rad, 0_rad)}};

  auto aStore = detail::makeAlignmentStore(
      gctx, detector.get(), globalShiftMean, globalAnglesStdErr,
      localShiftsMean, localShiftsStdErr, globalAnglesMean, globalAnglesStdErr,
      localAnglesMean, localAnglesStdErr);
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
  seqCfg.numThreads = 1;
  seqCfg.skip = 0;
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
  dummyReaderCfg.outputSourceLinkIndices = "SimMeasurementIndices";
  dummyReaderCfg.nEvents = 1e5;
  // dummyReaderCfg.nEvents = 3e3;

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

  // ClusterSizeBasedDigitizer::Config digitizerCfg;
  // digitizerCfg.clSizeProbsStdDevs = {{1, {0.168115, 0.00333208, 0.00395262}},
  //                                    {2, {0.284536, 0.00473498, 0.00444654}},
  //                                    {3, {0.217946, 0.00431227, 0.00417356}},
  //                                    {4, {0.329402, 0.00490848,
  //                                    0.00509384}}};
  // auto digitizer = std::make_shared<ClusterSizeBasedDigitizer>(digitizerCfg);

  // Vertex generator
  GaussianVertexGenerator::Config vertexGenCfg;
  vertexGenCfg.mean = 
    // Acts::Vector3(goInst.bpm0CenterPrimary - 5_mm, 0, 1_mm);
    Acts::Vector3(goInst.ipSurfaceCenterPrimary - 0.1_mm, 0, 1_mm);
  vertexGenCfg.cov = Acts::SquareMatrix3::Identity() * 30_um * 30_um;
  auto vertexGen = std::make_shared<GaussianVertexGenerator>(vertexGenCfg);

  SphericalMomentumGenerator::Config momGenCfg;
  momGenCfg.pRange = {3.0_GeV, 5.0_GeV};
  // momGenCfg.pRange = {2.0_GeV, 3.0_GeV};
  // momGenCfg.pRange = {1.9_GeV, 2.1_GeV};
  momGenCfg.phiRange = {-1e-3_rad, 1e-3_rad};
  momGenCfg.thetaRange = {M_PI_2 - 1e-3_rad, M_PI_2 + 1e-3_rad};

  auto momGen = std::make_shared<SphericalMomentumGenerator>(momGenCfg);

  // Measurement creator
  MeasurementsCreator::Config measCreatorCfg{
      .vertexGenerator = vertexGen,
      .momentumGenerator = momGen,
      .hitDigitizer = digitizer,
      .referenceSurface = refSurface.get(),
      .maxSteps = 1000,
      .isSignal = true,
      .hypothesis = Acts::ParticleHypothesis::electron(),
      .charge = -1_e};

  std::unordered_map<Acts::GeometryIdentifier, MeasurementsCreator::Constraints>
      measCreatorConstraints;
  for (auto& det : detector->detectorElements()) {
    const auto& surface = det->surface();
    const auto& geoId = surface.geometryId();
    if (geoId.sensitive() >= goInst.bpm0Parameters.geoId &&
        geoId.sensitive() <= goInst.bpm3Parameters.geoId) {
      measCreatorConstraints.insert({geoId, {-3e4, 3e4, -3e4, 3e4}});
    }
  }
  measCreatorCfg.constraints = measCreatorConstraints;

  auto measCreator = std::make_shared<MeasurementsCreator>(
      measCreatorPropagator, measCreatorCfg);

  MeasurementsEmbeddingAlgorithm::Config measCreatorAlgoCfg;
  measCreatorAlgoCfg.inputSourceLinks = "SimMeasurements";
  measCreatorAlgoCfg.inputSimClusters = "SimClusters";
  measCreatorAlgoCfg.inputSourceLinkIndices = "SimMeasurementIndices";
  measCreatorAlgoCfg.outputSourceLinks = "Measurements";
  measCreatorAlgoCfg.outputSimClusters = "Clusters";
  measCreatorAlgoCfg.outputSourceLinkIndices = "MeasurementIndices";
  measCreatorAlgoCfg.measurementGenerator = measCreator;
  measCreatorAlgoCfg.randomNumberSvc =
      std::make_shared<RandomNumbers>(RandomNumbers::Config());
  measCreatorAlgoCfg.nMeasurements = 1;

  sequencer.addAlgorithm(std::make_shared<MeasurementsEmbeddingAlgorithm>(
      measCreatorAlgoCfg, logLevel));

  // --------------------------------------------------------------
  // Add background

  // Background creator
  UniformBackgroundCreator::Config bkgCreatorCfg;
  bkgCreatorCfg.resolution = {5_um, 5_um};
  bkgCreatorCfg.nMeasurements = 700;
  bkgCreatorCfg.surfaces = detSurfaces;

  auto bkgCreator = std::make_shared<UniformBackgroundCreator>(bkgCreatorCfg);

  MeasurementsEmbeddingAlgorithm::Config bkgCreatorAlgoCfg;
  bkgCreatorAlgoCfg.inputSourceLinks = "Measurements1";
  bkgCreatorAlgoCfg.inputSimClusters = "Clusters1";
  bkgCreatorAlgoCfg.inputSourceLinkIndices = "MeasurementIndices1";
  bkgCreatorAlgoCfg.outputSourceLinks = "Measurements";
  bkgCreatorAlgoCfg.outputSimClusters = "Clusters";
  bkgCreatorAlgoCfg.outputSourceLinkIndices = "MeasurementIndices";
  bkgCreatorAlgoCfg.measurementGenerator = bkgCreator;
  bkgCreatorAlgoCfg.randomNumberSvc =
      std::make_shared<RandomNumbers>(RandomNumbers::Config());
  bkgCreatorAlgoCfg.nMeasurements = 1;

  // sequencer.addAlgorithm(std::make_shared<MeasurementsEmbeddingAlgorithm>(
  //     bkgCreatorAlgoCfg, logLevel));

  // --------------------------------------------------------------
  // Event write out

  // Sim cluster writer
  auto clusterWriterCfgSig = RootSimClusterWriter::Config();
  clusterWriterCfgSig.inputClusters = "Clusters";
  clusterWriterCfgSig.treeName = "clusters";
  clusterWriterCfgSig.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/sim/"
      "clusters.root";

  sequencer.addWriter(
      std::make_shared<RootSimClusterWriter>(clusterWriterCfgSig, logLevel));

  // Magnetic field writer
  auto magFieldWriterCfg = E320::E320MagneticFieldWriter::Config();
  magFieldWriterCfg.treeName = "magnets";
  magFieldWriterCfg.filePath =
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/sim/"
      "magnets.root";

  sequencer.addWriter(std::make_shared<E320::E320MagneticFieldWriter>(
      magFieldWriterCfg, logLevel));

  return sequencer.run();
}
