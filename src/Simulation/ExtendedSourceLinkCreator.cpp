#include "TrackingPipeline/Simulation/ExtendedSourceLinkCreator.hpp"

#include "TrackingPipeline/EventData/ExtendedSourceLink.hpp"

Acts::SourceLink ExtendedSourceLinkCreator::operator()(
    const Acts::GeometryContext& gctx, std::size_t eventNumber, std::size_t idx,
    const Acts::BoundTrackParameters& trackParameters,
    RandomEngine& rng) const {
  const Acts::BoundVector& boundVec = trackParameters.parameters();
  const Acts::Surface& refSurface = trackParameters.referenceSurface();
  const Acts::GeometryIdentifier& geoId = refSurface.geometryId();

  // Digitize hits
  Acts::Vector2 trueLocalPos(boundVec[Acts::eBoundLoc0],
                             boundVec[Acts::eBoundLoc1]);
  auto [digLocalPos, digHitCov] =
      m_cfg.hitDigitizer->digitize(rng, geoId, trueLocalPos);
  Acts::Vector3 digGlobalPos = trackParameters.referenceSurface().localToGlobal(
      gctx, digLocalPos, Acts::Vector3::UnitX());

  // Digitize directions
  Acts::Vector2 trueAngles(boundVec[Acts::eBoundPhi],
                           boundVec[Acts::eBoundTheta]);
  auto [digAngles, digAnglesCov] =
      m_cfg.angleDigitizer->digitize(rng, geoId, trueAngles);

  // Direction is specified in the track frame of reference,
  // which concides with the global (in angles)
  double phi = digAngles(0);
  double theta = digAngles(1);
  Acts::Vector3 digGlobalDir(std::sin(theta) * std::cos(phi),
                             std::sin(theta) * std::sin(phi), std::cos(theta));

  // Create measurement
  Acts::ActsVector<ExtendedSourceLink::localSubspaceSize> measLoc;
  measLoc.head(2) = digLocalPos;
  measLoc.tail(2) = digAngles;

  Acts::ActsVector<ExtendedSourceLink::globalSubspaceSize> measGlob;
  measGlob.head(3) = digGlobalPos;
  measGlob.tail(3) = digGlobalDir;

  Acts::ActsSquareMatrix<ExtendedSourceLink::localSubspaceSize> measCov =
      Acts::ActsSquareMatrix<ExtendedSourceLink::localSubspaceSize>::Zero();
  measCov.topLeftCorner(2, 2) = digHitCov;
  measCov.bottomRightCorner(2, 2) = digAnglesCov;

  ExtendedSourceLink esl(measLoc, measGlob, measCov, geoId, eventNumber, idx,
                         false);

  return Acts::SourceLink(esl);
}
