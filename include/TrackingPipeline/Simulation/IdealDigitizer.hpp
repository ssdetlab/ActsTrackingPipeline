#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

class IdealDigitizer : public IDigitizer {
 public:
  IdealDigitizer();

  Acts::Vector2 genCluster(RandomEngine& rng,
                           const Acts::GeometryIdentifier& geoId,
                           const Acts::Vector2& pos) const override;
  Acts::SquareMatrix2 getCovariance(
      const Acts::GeometryIdentifier& geoId) const override;

 private:
  Acts::SquareMatrix2 m_cov;
};
