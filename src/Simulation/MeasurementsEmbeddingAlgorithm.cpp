#include "TrackingPipeline/Simulation/MeasurementsEmbeddingAlgorithm.hpp"

#include <stdexcept>
#include <utility>

#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmRegistry.hpp"

#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"
#include "TrackingPipeline/Simulation/SimpleDigitizer.hpp"
#include "TrackingPipeline/Simulation/GaussianVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/SphericalMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/UniformBackgroundCreator.hpp"

#include "TrackingPipeline/TrackFitting/FittingServices.hpp"

#include <toml.hpp>


MeasurementsEmbeddingAlgorithm::MeasurementsEmbeddingAlgorithm(
    const Config& config, Acts::Logging::Level level)
    : IAlgorithm("MeasurementsEmbeddingAlgorithm", level), m_cfg(config) {
  if (m_cfg.measurementGenerator == nullptr) {
    throw std::runtime_error("Generator is not initialized");
  }

  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_inputSimClusters.initialize(m_cfg.inputSimClusters);

  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
  m_outputSimClusters.initialize(m_cfg.outputSimClusters);
}

ProcessCode MeasurementsEmbeddingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Get the inputs from the context
  auto sourceLinks = m_inputSourceLinks(ctx);
  auto clusters = m_inputSimClusters(ctx);

  ACTS_DEBUG("Received " << clusters.size() << " clusters");
  ACTS_DEBUG("Received " << sourceLinks.size() << " source links");

  RandomEngine rng = m_cfg.randomNumberSvc->spawnGenerator();

  // Create the measurements
  ACTS_DEBUG("Starting propagation of " << m_cfg.nMeasurements << " tracks");
  std::size_t inputSize = sourceLinks.size();
  std::size_t outputSize = inputSize + m_cfg.nMeasurements;
  sourceLinks.reserve(outputSize);
  int index = inputSize;
  for (std::size_t i = inputSize; i < outputSize; i++) {
    const auto& [sls, cls] = m_cfg.measurementGenerator->gen(ctx, rng, index);
    index += sls.size();

    ACTS_VERBOSE("Created " << sls.size() << " measurements");
    sourceLinks.insert(sourceLinks.end(), sls.begin(), sls.end());
    clusters.insert(clusters.end(), cls.begin(), cls.end());
  }

  ACTS_DEBUG("Created " << sourceLinks.size() << " measurements in total");
  m_outputSourceLinks(ctx, std::move(sourceLinks));
  m_outputSimClusters(ctx, std::move(clusters));

  return ProcessCode::SUCCESS;
}

// Registrars for SignalMeasurementEmbedding and BkgMeasurementEmbedding
namespace {

using TrackingPipeline::AlgorithmPtr;
using TrackingPipeline::AlgorithmRegistry;

// Signal: uses MeasurementsCreator with helpers from FittingServices
AlgorithmPtr buildSignalMeasurementEmbedding(
    const toml::value& section,
    Acts::Logging::Level logLevel) {

  auto& svc = FittingServices::instance();

  if (!svc.propagator || !svc.simDigitizer ||
      !svc.simVertexGenerator || !svc.simMomentumGenerator) {
    throw std::runtime_error(
        "SignalMeasurementEmbedding: FittingServices FastSim helpers not initialized");
  }

  // MeasurementCreator config: shared helpers + per-instance maxSteps
  MeasurementsCreator::Config mcfg{
      /*vertexGenerator*/   svc.simVertexGenerator,
      /*momentumGenerator*/ svc.simMomentumGenerator,
      /*hitDigitizer*/      svc.simDigitizer,
      /*referenceSurface*/  svc.referenceSurface.get(),
      /*maxSteps*/          static_cast<std::size_t>(
                                toml::find_or<int>(section, "meas_creator_maxSteps", 1000)),
      /*isSignal*/          true,
      /*hypothesis*/        Acts::ParticleHypothesis::electron(),  
      /*charge*/            -1.0,                              
      /*constraints*/       {}                                 
  };

  auto measCreator =
      std::make_shared<MeasurementsCreator>(*svc.propagator, mcfg);

  // Algorithm config from TOML
  MeasurementsEmbeddingAlgorithm::Config cfg;
  cfg.inputSourceLinks  = toml::find<std::string>(section, "inputSourceLinks");
  cfg.inputSimClusters  = toml::find<std::string>(section, "inputSimClusters");
  cfg.outputSourceLinks = toml::find<std::string>(section, "outputSourceLinks");
  cfg.outputSimClusters = toml::find<std::string>(section, "outputSimClusters");
  cfg.nMeasurements     = toml::find_or<std::size_t>(section, "nMeasurements", 0);

  cfg.measurementGenerator = measCreator;
  cfg.randomNumberSvc      =
      std::make_shared<RandomNumbers>(RandomNumbers::Config());

  return std::make_shared<MeasurementsEmbeddingAlgorithm>(cfg, logLevel);
}

// Background: uses UniformBackgroundCreator with surfaces from gDetSurfaces
AlgorithmPtr buildBackgroundMeasurementEmbedding(
    const toml::value& section,
    Acts::Logging::Level logLevel) {

  auto& svc = FittingServices::instance();
  if (svc.simDetSurfaces.empty()) {
    throw std::runtime_error(
        "BkgMeasurementEmbedding: simDetSurfaces is empty; detector surfaces were not initialized");
  }

  UniformBackgroundCreator::Config bcfg;

  // Resolution as [dx, dy]; default 5e-6 if not provided
  auto res = toml::find_or<std::vector<double>>(
      section, "bkgCreator_resolution", {5.0e-6, 5.0e-6});
  if (res.size() != 2u) {
    throw std::runtime_error(
        "BkgMeasurementEmbedding: bkgCreator_resolution must be [dx, dy]");
  }
  bcfg.resolution    = {res[0], res[1]};
  bcfg.nMeasurements =
      toml::find_or<int>(section, "bkgCreator_nMeasurements", 20);
  bcfg.surfaces      = svc.simDetSurfaces;

  auto bkgCreator = std::make_shared<UniformBackgroundCreator>(bcfg);

  MeasurementsEmbeddingAlgorithm::Config cfg;
  cfg.inputSourceLinks  = toml::find<std::string>(section, "inputSourceLinks");
  cfg.inputSimClusters  = toml::find<std::string>(section, "inputSimClusters");
  cfg.outputSourceLinks = toml::find<std::string>(section, "outputSourceLinks");
  cfg.outputSimClusters = toml::find<std::string>(section, "outputSimClusters");
  cfg.nMeasurements     = toml::find_or<std::size_t>(section, "nMeasurements", 0);

  cfg.measurementGenerator = bkgCreator;
  cfg.randomNumberSvc      =
      std::make_shared<RandomNumbers>(RandomNumbers::Config());

  return std::make_shared<MeasurementsEmbeddingAlgorithm>(cfg, logLevel);
}

struct MeasurementsEmbeddingRegistrar {
  MeasurementsEmbeddingRegistrar() {
    // Factory keys used from TOML "type" fields
    AlgorithmRegistry::instance().registerBuilder(
        "SignalMeasurementEmbedding", buildSignalMeasurementEmbedding);
    AlgorithmRegistry::instance().registerBuilder(
        "BkgMeasurementEmbedding", buildBackgroundMeasurementEmbedding);
  }
};

MeasurementsEmbeddingRegistrar _MeasurementsEmbeddingRegistrar;

}  // namespace