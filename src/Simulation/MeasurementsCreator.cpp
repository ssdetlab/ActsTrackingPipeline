#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/VectorHelpers.hpp>

#include <cstddef>
#include <cstdlib>
#include <stdexcept>
#include <utility>
#include <vector>

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreatorAction.hpp"

MeasurementsCreator::MeasurementsCreator(const Propagator propagator,
                                         const Config& config)
    : m_cfg(config), m_propagator(propagator) {
  Acts::BoundVector ipStdDev;
  ipStdDev[Acts::eBoundLoc0] = 100_um;
  ipStdDev[Acts::eBoundLoc1] = 100_um;
  ipStdDev[Acts::eBoundTime] = 25_ns;
  ipStdDev[Acts::eBoundPhi] = 2_degree;
  ipStdDev[Acts::eBoundTheta] = 2_degree;
  ipStdDev[Acts::eBoundQOverP] = 1 / 100_GeV;
  m_ipCov = ipStdDev.cwiseProduct(ipStdDev).asDiagonal();
};

std::tuple<std::vector<Acts::SourceLink>, SimClusters> MeasurementsCreator::gen(
    const AlgorithmContext& ctx, RandomEngine& rng, std::size_t id) const {
  using Actions = Acts::ActionList<MeasurementsCreatorAction>;
  using Aborters = Acts::AbortList<Acts::EndOfWorldReached>;
  using PropagatorOptions =
      typename Propagator::template Options<Actions, Aborters>;

  // Set options for propagator
  PropagatorOptions options(ctx.geoContext, ctx.magFieldContext);

  auto& creator = options.actionList.template get<MeasurementsCreatorAction>();
  options.maxSteps = m_cfg.maxSteps;
  creator.sourceId = id;

  int index = static_cast<int>(id);

  Acts::Vector3 spatial = m_cfg.vertexGenerator->genVertex(rng);
  Acts::Vector4 mPos4 = {spatial.x(), spatial.y(), spatial.z(), 0};

  Acts::Vector3 mom = m_cfg.momentumGenerator->genMomentum(rng);
  double p = mom.norm();
  double phi = Acts::VectorHelpers::phi(mom);
  double theta = Acts::VectorHelpers::theta(mom);

  TrackParameters trackParameters(mPos4, phi, theta, m_cfg.charge / p, m_ipCov,
                                  m_cfg.hypothesis);

  MeasurementsCreatorAction::result_type resultParameters;
  try {
    auto result =
        m_propagator.propagate(std::move(trackParameters), std::move(options))
            .value();

    resultParameters =
        result.template get<MeasurementsCreatorAction::result_type>();
  } catch (const std::runtime_error& err) {
    std::cout << err.what() << "\n";
  }

  int trackId = (m_cfg.isSignal) ? 1 : -1;
  std::size_t resSize = resultParameters.size();

  SimClusters simClusters;
  simClusters.reserve(resSize);
  std::vector<Acts::SourceLink> sourceLinks;
  sourceLinks.reserve(resSize);

  for (std::size_t i = 0; i < resSize; i++) {
    const auto& boundPars = resultParameters.at(i);
    const Acts::BoundVector& boundVec = boundPars.parameters();

    Acts::GeometryIdentifier geoId = boundPars.referenceSurface().geometryId();

    // Sometimes scattering makes particles
    // to be stuck in the same surface
    for (int j = simClusters.size() - 1; j >= 0; j--) {
      if (simClusters.at(j).sourceLink.geometryId() == geoId) {
        continue;
      }
    }

    // Digitize hits
    Acts::Vector2 trueLocalPos{boundVec[Acts::eBoundLoc0],
                               boundVec[Acts::eBoundLoc1]};
    Acts::Vector3 trueGlobalPos = boundPars.referenceSurface().localToGlobal(
        ctx.geoContext, trueLocalPos, Acts::Vector3::UnitX());

    auto [digCov, digLocalPos] =
        m_cfg.hitDigitizer->genCluster(rng, geoId, trueLocalPos);
    Acts::Vector3 digGlobalPos = boundPars.referenceSurface().localToGlobal(
        ctx.geoContext, digLocalPos, Acts::Vector3::UnitX());

    // Truth information
    Acts::SourceLink hitSl(SimpleSourceLink(trueLocalPos, trueGlobalPos, digCov,
                                            geoId, ctx.eventNumber, index + i));

    SimHit sm{
        boundVec, trueGlobalPos,        trackParameters,
        trackId,  static_cast<int>(id), static_cast<int>(ctx.eventNumber)};

    // Observable information
    SimpleSourceLink simpleSl(digLocalPos, digGlobalPos, digCov, geoId,
                              ctx.eventNumber, index + i);
    sourceLinks.push_back(Acts::SourceLink(simpleSl));

    SimCluster cl{
        simpleSl,
        {sm},
        m_cfg.isSignal,
    };
    simClusters.push_back(cl);
  }

  return {std::move(sourceLinks), std::move(simClusters)};
};
