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
#include <vector>

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"
#include "TrackingPipeline/Simulation/MeasurementsCreatorAction.hpp"

MeasurementsCreator::MeasurementsCreator(const Propagator& propagator,
                                         const Config& cfg)
    : m_cfg(cfg), m_propagator(propagator) {
  m_freeIpCov = Acts::FreeMatrix::Zero();

  m_freeIpCov.block(Acts::eFreePos0, Acts::eFreePos0, 3, 3) =
      m_cfg.vertexGenerator->getVertexCovariance();

  m_freeIpCov(Acts::eFreeTime, Acts::eFreeTime) = 25_ns;

  m_freeIpCov.block(Acts::eFreeDir0, Acts::eFreeDir0, 4, 4) =
      m_cfg.momentumGenerator->getMomentumCovariance();
};

std::size_t MeasurementsCreator::gen(
    const AlgorithmContext& ctx, RandomEngine& rng, std::size_t id,
    std::vector<Acts::SourceLink>& sourceLinks, SimClusters& simClusters,
    std::vector<std::size_t>& sourceLinksIndices) const {
  using Actions = Acts::ActionList<MeasurementsCreatorAction>;
  using Aborters = Acts::AbortList<Acts::EndOfWorldReached>;
  using PropagatorOptions =
      typename Propagator::template Options<Actions, Aborters>;

  // Set options for propagator
  PropagatorOptions options(ctx.geoContext, ctx.magFieldContext);

  auto& creator = options.actionList.template get<MeasurementsCreatorAction>();
  options.maxSteps = m_cfg.maxSteps;
  creator.sourceId = id;

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

  MeasurementsCreatorAction::result_type resultParameters;
  try {
    auto result = m_propagator.propagate(trackParameters, options).value();

    resultParameters =
        result.template get<MeasurementsCreatorAction::result_type>();
  } catch (const std::runtime_error& err) {
    std::cout << err.what() << "\n";
  }

  int trackId = (m_cfg.isSignal) ? 1 : -1;
  std::size_t currentSize = sourceLinksIndices.size();
  std::size_t resSize = resultParameters.size();

  simClusters.reserve(currentSize + resSize);
  sourceLinks.reserve(currentSize + resSize);
  sourceLinksIndices.reserve(currentSize + resSize);

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

    // Observable data
    Acts::SourceLink sl = m_cfg.sourceLinkCreators.at(geoId.sensitive())(
        ctx.geoContext, ctx.eventNumber, sourceLinks.size(), boundPars, rng);
    sourceLinksIndices.push_back(sourceLinks.size());
    sourceLinks.push_back(sl);

    // Truth information
    Acts::Vector2 trueLocalPos{boundVec[Acts::eBoundLoc0],
                               boundVec[Acts::eBoundLoc1]};
    Acts::Vector3 trueGlobalPos = boundPars.referenceSurface().localToGlobal(
        ctx.geoContext, trueLocalPos, Acts::Vector3::UnitX());
    SimHit sh{
        boundVec, trueGlobalPos,        trackParameters,
        trackId,  static_cast<int>(id), static_cast<int>(ctx.eventNumber)};
    SimCluster cl{
        sl,
        {sh},
        m_cfg.isSignal,
    };
    simClusters.push_back(cl);
  }

  return resSize;
};
