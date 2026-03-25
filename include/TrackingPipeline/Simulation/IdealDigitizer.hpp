#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

class IdealDigitizer : public IDigitizer {
 public:
  IdealDigitizer();

  std::pair<Acts::Vector2, Acts::SquareMatrix2> genCluster(
      RandomEngine& rng, const Acts::GeometryIdentifier& geoId,
      const Acts::Vector2& pos) const override;

 private:
  Acts::SquareMatrix2 m_cov;
};
