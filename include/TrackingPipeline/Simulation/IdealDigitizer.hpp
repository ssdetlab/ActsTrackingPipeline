#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

struct IdealDigitizer : public IDigitizer {
  std::pair<Acts::SquareMatrix2, Acts::Vector2> genCluster(
      RandomEngine& rng, Acts::GeometryIdentifier /*geoId*/,
      Acts::Vector2 pos) const override {
    Acts::SquareMatrix2 cov = Acts::SquareMatrix2::Zero();

    return {cov, pos};
  }
};
