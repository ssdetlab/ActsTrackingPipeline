#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
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
#include "TrackingPipeline/Io/DummyReader.hpp"
#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"
#include "TrackingPipeline/Simulation/GaussianVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"
#include "TrackingPipeline/Simulation/SimpleDigitizer.hpp"
#include "TrackingPipeline/Simulation/SphericalMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/UniformBackgroundCreator.hpp"

using ActionList = Acts::ActionList<>;
using AbortList = Acts::AbortList<Acts::EndOfWorldReached>;

using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                    Acts::Experimental::DetectorNavigator>;
using PropagatorOptions =
    typename Propagator::template Options<ActionList, AbortList>;

using RecoTrajectory = Acts::VectorMultiTrajectory;
using RecoTrackContainer = Acts::VectorTrackContainer;
using KF = Acts::KalmanFitter<Propagator, RecoTrajectory>;

using CKFTrackContainer = Acts::TrackContainer<Acts::VectorTrackContainer,
                                               Acts::VectorMultiTrajectory,
                                               Acts::detail::ValueHolder>;

using TrackStateContainerBackend =
    typename CKFTrackContainer::TrackStateContainerBackend;

using namespace Acts::UnitLiterals;

namespace eg = E320Geometry;

std::unique_ptr<const eg::GeometryOptions> eg::GeometryOptions::m_instance =
    nullptr;

int main(int argc, char* argv[]) {
  int id = std::stoi(argv[1]);

  const auto& goInst = *eg::GeometryOptions::instance();

  // Set the log level
  Acts::Logging::Level logLevel = Acts::Logging::FATAL;

  // Dummy context and options
  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext cctx;
  E320Geometry::GeometryOptions gOpt;

  // --------------------------------------------------------------
  // Detector setup

  auto detector = eg::buildDetector(gctx);

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
      if (surf->geometryId().sensitive()) {
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
      {10, Acts::Vector3(0_mm, 30_um, 30_um)},
      {12, Acts::Vector3(0_mm, 30_um, 30_um)},
      {14, Acts::Vector3(0_mm, 30_um, 30_um)},
      {16, Acts::Vector3(0_mm, 30_um, 30_um)},
      {18, Acts::Vector3(0_mm, 30_um, 30_um)}};

  Acts::Vector3 globalAnglesMean(0_rad, 0_rad, 0_rad);
  Acts::Vector3 globalAnglesStdErr(0_rad, 0_rad, 0_rad);

  std::unordered_map<int, Acts::Vector3> localAnglesMean{
      {10, Acts::Vector3(0_rad, 0_rad, 0_rad)},
      {12, Acts::Vector3(0_rad, 0_rad, 0_rad)},
      {14, Acts::Vector3(0_rad, 0_rad, 0_rad)},
      {16, Acts::Vector3(0_rad, 0_rad, 0_rad)},
      {18, Acts::Vector3(0_rad, 0_rad, 0_rad)}};
  std::unordered_map<int, Acts::Vector3> localAnglesStdErr{
      {10, Acts::Vector3(0_rad, 0_rad, 1e-3_rad)},
      {12, Acts::Vector3(0_rad, 0_rad, 1e-3_rad)},
      {14, Acts::Vector3(0_rad, 0_rad, 1e-3_rad)},
      {16, Acts::Vector3(0_rad, 0_rad, 1e-3_rad)},
      {18, Acts::Vector3(0_rad, 0_rad, 1e-3_rad)}};

  auto aStore = detail::makeAlignmentStore(
      gctx, detector.get(), globalShiftMean, globalAnglesStdErr,
      localShiftsMean, localShiftsStdErr, globalAnglesMean, globalAnglesStdErr,
      localAnglesMean, localAnglesStdErr);
  AlignmentContext alignCtx(aStore);
  Acts::GeometryContext testCtx{alignCtx};
  for (auto& v : detector->volumes()) {
    for (auto& s : v->surfaces()) {
      if (s->geometryId().sensitive()) {
        std::cout << "-----------------------------------\n";
        std::cout << "SURFACE " << s->geometryId() << "\n";
        std::cout << "CENTER " << s->center(testCtx).transpose() << " -- "
                  << s->center(Acts::GeometryContext()).transpose() << "\n";
        std::cout << "NORMAL "
                  << s->normal(testCtx, s->center(testCtx),
                               Acts::Vector3::UnitY())
                         .transpose()
                  << " -- "
                  << s->normal(Acts::GeometryContext(),
                               s->center(Acts::GeometryContext()),
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

  auto field = eg::buildMagField(gctx);

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
  refSurface->assignGeometryId(std::move(geoId));

  // --------------------------------------------------------------
  // Event reading
  SimpleSourceLink::SurfaceAccessor surfaceAccessor{detector.get()};

  // Setup the sequencer
  Sequencer::Config seqCfg;
  seqCfg.numThreads = 1;
  seqCfg.skip = 100000 * id;
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
  dummyReaderCfg.nEvents = (id + 1) * 100000;

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
  GaussianVertexGenerator::Config vertexGenCfg;
  vertexGenCfg.mean = Acts::Vector3(0, 0, 0);
  vertexGenCfg.cov = Acts::SquareMatrix3::Identity() * 30_um;
  auto vertexGen = std::make_shared<GaussianVertexGenerator>(vertexGenCfg);

  SphericalMomentumGenerator::Config momGenCfg;
  momGenCfg.pRange = {2.0_GeV, 3.0_GeV};
  momGenCfg.phiRange = {-0.001, 0.001};
  momGenCfg.thetaRange = {M_PI_2 - 0.003, M_PI_2 + 0.003};

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
    if (geoId.sensitive() && geoId.sensitive() >= 40) {
      measCreatorConstraints.insert({geoId, {-3, 3, -3, 3}});
    }
  }
  measCreatorCfg.constraints = measCreatorConstraints;

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
      "clusters-" +
      std::to_string(id) + ".root";

  sequencer.addWriter(
      std::make_shared<RootSimClusterWriter>(clusterWriterCfgSig, logLevel));

  return sequencer.run();
}
