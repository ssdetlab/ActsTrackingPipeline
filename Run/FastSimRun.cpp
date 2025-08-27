#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <Eigen/src/Core/ArithmeticSequence.h>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/ApollonGeometry.hpp"
#include "TrackingPipeline/Geometry/ApollonGeometryConstraints.hpp"
#include "TrackingPipeline/Infrastructure/Sequencer.hpp"
#include "TrackingPipeline/Io/DummyReader.hpp"
#include "TrackingPipeline/Io/RootSimClusterWriter.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"
#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"
#include "TrackingPipeline/Simulation/SimpleDigitizer.hpp"
#include "TrackingPipeline/Simulation/SphericalMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/StationaryMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/StationaryVertexGenerator.hpp"
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

namespace ag = ApollonGeometry;

std::unique_ptr<const ag::GeometryOptions> ag::GeometryOptions::m_instance =
    nullptr;

int main() {
  // Set the log level
  Acts::Logging::Level logLevel = Acts::Logging::INFO;

  // Context and options
  Acts::GeometryContext gctx;

  // --------------------------------------------------------------
  // Detector setup

  auto detector = ApollonGeometry::buildDetector(gctx);

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

  auto field = ApollonGeometry::buildMagField(gctx);

  // --------------------------------------------------------------
  // Event reading
  SimpleSourceLink::SurfaceAccessor surfaceAccessor{detector.get()};

  // Setup the sequencer
  Sequencer::Config seqCfg;
  seqCfg.numThreads = 1;
  seqCfg.skip = 1;
  seqCfg.trackFpes = false;
  seqCfg.logLevel = logLevel;
  Sequencer sequencer(seqCfg);

  // --------------------------------------------------------------
  // Add dummy reader
  DummyReader::Config dummyReaderCfg;
  dummyReaderCfg.outputSourceLinks = "SimMeasurements";
  dummyReaderCfg.outputSimClusters = "SimClusters";
  dummyReaderCfg.nEvents = 1e3 + 1;

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
  auto vertexGen = std::make_shared<StationaryVertexGenerator>();
  vertexGen->vertex = Acts::Vector3(-2_mm, 0, 0);

  // Momentum generator
  // auto momGen = std::make_shared<SphericalMomentumGenerator>();
  // momGen->pRange = {0.3_GeV, 0.7_GeV};
  // momGen->phiRange = {0, 2 * M_PI};
  // momGen->thetaRange = {M_PI_2 - 3e-3, M_PI_2 + 3e-3};

  auto momGen = std::make_shared<StationaryMomentumGenerator>();
  momGen->momentum = 0.7_GeV * Acts::Vector3::UnitX();

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
  measCreatorAlgoCfg.outputSourceLinks = "Measurements1";
  measCreatorAlgoCfg.outputSimClusters = "Clusters1";
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
  bkgCreatorCfg.nMeasurements = 20;
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

  sequencer.addAlgorithm(std::make_shared<MeasurementsEmbeddingAlgorithm>(
      bkgCreatorAlgoCfg, logLevel));

  // --------------------------------------------------------------
  // Event write out

  // Sim cluster writer
  auto clusterWriterCfgSig = RootSimClusterWriter::Config();
  clusterWriterCfgSig.inputClusters = "Clusters1";
  clusterWriterCfgSig.treeName = "clusters";
  clusterWriterCfgSig.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/fast_sim_data/"
      "clusters-sim-sig.root";

  sequencer.addWriter(
      std::make_shared<RootSimClusterWriter>(clusterWriterCfgSig, logLevel));

  auto clusterWriterCfgSigBkg = RootSimClusterWriter::Config();
  clusterWriterCfgSigBkg.inputClusters = "Clusters";
  clusterWriterCfgSigBkg.treeName = "clusters";
  clusterWriterCfgSigBkg.filePath =
      "/home/romanurmanov/work/Apollon/tracking/out_data/fast_sim_data/"
      "clusters-sim-sig-bkg.root";

  sequencer.addWriter(
      std::make_shared<RootSimClusterWriter>(clusterWriterCfgSig, logLevel));

  sequencer.addWriter(
      std::make_shared<RootSimClusterWriter>(clusterWriterCfgSigBkg, logLevel));

  return sequencer.run();
}
