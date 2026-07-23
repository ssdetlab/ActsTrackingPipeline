#include "TrackingPipeline/Simulation/SimpleSourceLinkCreator.hpp"

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

Acts::SourceLink SimpleSourceLinkCreator::operator()(
    const Acts::GeometryContext& gctx, std::size_t eventNumber, std::size_t idx,
    const Acts::BoundTrackParameters& trackParameters,
    RandomEngine& rng) const {
  const Acts::BoundVector& boundVec = trackParameters.parameters();
  const Acts::Surface& refSurface = trackParameters.referenceSurface();
  const Acts::GeometryIdentifier& geoId = refSurface.geometryId();

  // Digitize hits
  Acts::Vector2 trueLocalPos(boundVec[Acts::eBoundLoc0],
                             boundVec[Acts::eBoundLoc1]);
  auto [digLocalPos, digCov] =
      m_cfg.hitDigitizer->digitize(rng, geoId, trueLocalPos);
  Acts::Vector3 digGlobalPos = trackParameters.referenceSurface().localToGlobal(
      gctx, digLocalPos, Acts::Vector3::UnitX());

  // Create measurement
  SimpleSourceLink ssl(digLocalPos, digGlobalPos, digCov, geoId, eventNumber,
                       idx);

  return Acts::SourceLink(ssl);
}
