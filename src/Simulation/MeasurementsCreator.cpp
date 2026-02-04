#include "TrackingPipeline/Simulation/MeasurementsCreator.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreatorAction.hpp"

MeasurementsCreator::MeasurementsCreator(const Propagator propagator,
                                         const Config& cfg)
    : m_cfg(cfg), m_propagator(propagator) {
  m_freeIpCov = Acts::FreeMatrix::Zero();

  m_freeIpCov.block(Acts::eFreePos0, Acts::eFreePos0, 3, 3) =
      m_cfg.vertexGenerator->getCovariance();

  m_freeIpCov(Acts::eFreeTime, Acts::eFreeTime) = 25_ns;

  m_freeIpCov.block(Acts::eFreeDir0, Acts::eFreeDir0, 4, 4) =
      m_cfg.momentumGenerator->getCovariance();
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

  Acts::FreeToBoundMatrix jacToLoc =
      m_cfg.referenceSurface->freeToBoundJacobian(ctx.geoContext, spatial,
                                                  mom.normalized());
  Acts::BoundMatrix ipCov = jacToLoc * m_freeIpCov * jacToLoc.transpose();

  TrackParameters trackParameters(mPos4, phi, theta, m_cfg.charge / p, ipCov,
                                  m_cfg.hypothesis);
  trackParameters.absoluteMomentum();

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
    const Acts::BoundTrackParameters& boundPars = resultParameters.at(i);
    const Acts::BoundVector& boundVec = boundPars.parameters();

    Acts::GeometryIdentifier geoId = boundPars.referenceSurface().geometryId();
    if (m_cfg.constraints.contains(geoId)) {
      if (boundVec[Acts::eBoundLoc0] < m_cfg.constraints.at(geoId).minLocX ||
          boundVec[Acts::eBoundLoc0] > m_cfg.constraints.at(geoId).maxLocX ||
          boundVec[Acts::eBoundLoc1] < m_cfg.constraints.at(geoId).minLocY ||
          boundVec[Acts::eBoundLoc1] > m_cfg.constraints.at(geoId).maxLocY) {
        return {};
      }
    }

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

    Acts::Vector2 digLocalPos =
        m_cfg.hitDigitizer->genCluster(rng, geoId, trueLocalPos);
    Acts::Vector3 digGlobalPos = boundPars.referenceSurface().localToGlobal(
        ctx.geoContext, digLocalPos, Acts::Vector3::UnitX());
    Acts::SquareMatrix2 digCov = m_cfg.hitDigitizer->getCovariance(geoId);

    // Truth information
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
