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
#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"
#include "TrackingPipeline/MagneticField/ConstantMagField.hpp"
#include "TrackingPipeline/MagneticField/IdealQuadrupoleMagField.hpp"
#include "TrackingPipeline/MagneticField/MagneticFieldContextDecorator.hpp"
#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"
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
  Acts::Vector3 globalShiftMean(0, -4_mm, 5_mm);
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

  Acts::Vector3 globalAnglesMean(0_rad, 0_rad, 2e-3_rad);
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

  auto mStore1 = std::make_shared<MagneticFieldStore>();
  mStore1->store = {
      {goInst.quad1Id,
       Acts::MagneticFieldProvider::Cache(
           std::in_place_type<IdealQuadrupoleMagField::Cache>, 0)},
      {goInst.quad2Id,
       Acts::MagneticFieldProvider::Cache(
           std::in_place_type<IdealQuadrupoleMagField::Cache>, 0)},
      {goInst.quad3Id,
       Acts::MagneticFieldProvider::Cache(
           std::in_place_type<IdealQuadrupoleMagField::Cache>, 0)},
      {goInst.xCorrectorId, Acts::MagneticFieldProvider::Cache(
                                std::in_place_type<ConstantMagField::Cache>,
                                Acts::Vector3(0, 0.026107_T, 0))},
      {goInst.dipoleId, Acts::MagneticFieldProvider::Cache(
                            std::in_place_type<ConstantMagField::Cache>,
                            Acts::Vector3(0, 0, -0.2192_T))}};

  auto mStore2 = std::make_shared<MagneticFieldStore>();
  mStore2->store = {
      {goInst.quad1Id,
       Acts::MagneticFieldProvider::Cache(
           std::in_place_type<IdealQuadrupoleMagField::Cache>, 0)},
      {goInst.quad2Id,
       Acts::MagneticFieldProvider::Cache(
           std::in_place_type<IdealQuadrupoleMagField::Cache>, 0)},
      {goInst.quad3Id,
       Acts::MagneticFieldProvider::Cache(
           std::in_place_type<IdealQuadrupoleMagField::Cache>, 0)},
      {goInst.xCorrectorId, Acts::MagneticFieldProvider::Cache(
                                std::in_place_type<ConstantMagField::Cache>,
                                Acts::Vector3(0, 0.026107_T, 0))},
      {goInst.dipoleId, Acts::MagneticFieldProvider::Cache(
                            std::in_place_type<ConstantMagField::Cache>,
                            Acts::Vector3(0, 0, -0.2192_T))}};

  MagneticFieldStoreCollection mFieldStoreCollection = {{0, mStore1},
                                                        {1000, mStore2}};
  auto mFieldParametersContext =
      std::make_shared<MagneticFieldParametersContext>(mFieldStoreCollection);

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
  sequencer.addContextDecorator(
      std::make_shared<MagneticFieldContextDecorator>(mFieldParametersContext));

  // --------------------------------------------------------------
  // Add dummy reader
  DummyReader::Config dummyReaderCfg;
  dummyReaderCfg.outputSourceLinks = "SimMeasurements";
  dummyReaderCfg.outputSimClusters = "SimClusters";
  dummyReaderCfg.outputSourceLinkIndices = "SimMeasurementIndices";
  dummyReaderCfg.nEvents = 2e3;

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

  SurfaceRangedDigitizer::Config digitizerCfg;
  for (const auto& surf : detSurfaces) {
    SurfaceRangedDigitizer::Resolution res =
        (geoId.sensitive() < goInst.bpm0Parameters.geoId)
            ? std::make_pair(5_um, 5_um)
            : std::make_pair(100_um, 100_um);
    digitizerCfg.resolutions.insert({surf->geometryId(), res});
  }
  auto digitizer = std::make_shared<SurfaceRangedDigitizer>(digitizerCfg);

  // ClusterSizeBasedDigitizer::Config digitizerCfg;
  // digitizerCfg.clSizeProbsStdDevs = {{1, {0.168115, 0.00333208, 0.00395262}},
  //                                    {2, {0.284536, 0.00473498, 0.00444654}},
  //                                    {3, {0.217946, 0.00431227, 0.00417356}},
  //                                    {4, {0.329402, 0.00490848,
  //                                    0.00509384}}};
  // auto digitizer = std::make_shared<ClusterSizeBasedDigitizer>(digitizerCfg);

  // Vertex generator
  GaussianVertexGenerator::Config vertexGenCfg;
  vertexGenCfg.mean = Acts::Vector3(goInst.bpm0CenterPrimary - 1_mm, 0, 0);
  vertexGenCfg.cov = Acts::SquareMatrix3::Identity() * 0_um;
  auto vertexGen = std::make_shared<GaussianVertexGenerator>(vertexGenCfg);

  SphericalMomentumGenerator::Config momGenCfg;
  momGenCfg.pRange = {2.5_GeV, 2.5_GeV};
  momGenCfg.phiRange = {0.0, 0.0};
  momGenCfg.thetaRange = {M_PI_2 + 0.0, M_PI_2 + 0.0};

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
      .charge = 1_e};

  std::unordered_map<Acts::GeometryIdentifier, MeasurementsCreator::Constraints>
      measCreatorConstraints;
  for (auto& det : detector->detectorElements()) {
    const auto& surface = det->surface();
    const auto& geoId = surface.geometryId();
    if (geoId.sensitive() != 0u &&
        geoId.sensitive() >= goInst.bpm0Parameters.geoId) {
      measCreatorConstraints.insert({geoId, {-30, 30, -30, 30}});
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
  bkgCreatorCfg.nMeasurements = 300;
  bkgCreatorCfg.surfaces = detSurfaces;

  auto bkgCreator = std::make_shared<UniformBackgroundCreator>(bkgCreatorCfg);

  MeasurementsEmbeddingAlgorithm::Config bkgCreatorAlgoCfg;
  bkgCreatorAlgoCfg.inputSourceLinks = "Measurements1";
  bkgCreatorAlgoCfg.inputSimClusters = "Clusters1";
  bkgCreatorAlgoCfg.outputSourceLinks = "Measurements";
  bkgCreatorAlgoCfg.outputSimClusters = "Clusters";
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

  return sequencer.run();
}
