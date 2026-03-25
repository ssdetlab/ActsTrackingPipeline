#include "TrackingPipeline/Simulation/IdealDigitizer.hpp"

IdealDigitizer::IdealDigitizer() {
  m_cov = Acts::SquareMatrix2::Zero();
}

std::pair<Acts::Vector2, Acts::SquareMatrix2> IdealDigitizer::genCluster(
    RandomEngine& /*rng*/, const Acts::GeometryIdentifier& /*geoId*/,
    const Acts::Vector2& pos) const {
  return {pos, m_cov};
}
