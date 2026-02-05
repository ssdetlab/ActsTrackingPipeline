#include "TrackingPipeline/Simulation/IdealDigitizer.hpp"

IdealDigitizer::IdealDigitizer() {
  m_cov = Acts::SquareMatrix2::Zero();
}

Acts::Vector2 IdealDigitizer::genCluster(
    RandomEngine& /*rng*/, const Acts::GeometryIdentifier& /*geoId*/,
    const Acts::Vector2& pos) const {
  return pos;
}

Acts::SquareMatrix2 IdealDigitizer::getCovariance(
    const Acts::GeometryIdentifier& /*geoId*/) const {
  return m_cov;
}
