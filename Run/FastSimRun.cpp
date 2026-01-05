#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <memory>
#include <vector>

#include "TrackingPipeline/Geometry/E320Geometry.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/DummyReader.hpp"
#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"
#include "TrackingPipeline/Simulation/GaussianVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"
#include "TrackingPipeline/Simulation/SimpleDigitizer.hpp"
#include "TrackingPipeline/Simulation/StationaryMomentumGenerator.hpp"

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

int main() {
  const auto& goInst = *eg::GeometryOptions::instance();

  // Set the log level
  Acts::Logging::Level logLevel = Acts::Logging::DEBUG;

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
  // The magnetic field setup

  auto field = eg::buildMagField(gctx);

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

  // --------------------------------------------------------------
  // Add dummy reader
  DummyReader::Config dummyReaderCfg;
  dummyReaderCfg.outputSourceLinks = "SimMeasurements";
  dummyReaderCfg.outputSimClusters = "SimClusters";
  dummyReaderCfg.nEvents = 10000;

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
  Acts::Vector3 vertexMean = Acts::Vector3::Zero();
  Acts::SquareMatrix3 vertexCov = Acts::SquareMatrix3::Identity() * 10_um;
  auto vertexGen =
      std::make_shared<GaussianVertexGenerator>(vertexMean, vertexCov);

  auto momGen = std::make_shared<StationaryMomentumGenerator>();
  momGen->momentum = {2.5_GeV, 0, 0};

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
  // Event write out

  // Sim cluster writer
  auto clusterWriterCfgSig = RootSimClusterWriter::Config();
  clusterWriterCfgSig.inputClusters = "Clusters";
  clusterWriterCfgSig.treeName = "clusters";
  clusterWriterCfgSig.filePath =
      "/home/romanurmanov/work/TrackingPipeline/ActsTrackingPipeline_build/"
      "clusters.root";

  sequencer.addWriter(
      std::make_shared<RootSimClusterWriter>(clusterWriterCfgSig, logLevel));

  return sequencer.run();
}
