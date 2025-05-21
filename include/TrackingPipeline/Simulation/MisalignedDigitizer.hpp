#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <map>

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

/// @brief Class that digitizes hits based on the provided
/// resolution assuming Gaussian smearing
struct MisalignedDigitizer : public IDigitizer {
  std::map<int, std::pair<double, double>> shifts;
  std::pair<double, double> resolution;

  std::pair<Acts::SquareMatrix2, Acts::Vector2> genCluster(
      RandomEngine& rng, Acts::GeometryIdentifier geoId,
      Acts::Vector2 pos) const override {
    std::normal_distribution<double> normalDist(0., 1.);

    int staveId = static_cast<int>(geoId.sensitive() - 1);
    Acts::Vector2 shift = {shifts.at(staveId).first, shifts.at(staveId).second};
    Acts::Vector2 res = {resolution.first, resolution.second};
    Acts::Vector2 digLocal =
        pos + shift +
        res.cwiseProduct(Acts::Vector2(normalDist(rng), normalDist(rng)));
    Acts::Vector2 stdDev = {
        /*std::hypot(resolution.first, 1e2 * shifts.at(staveId).first),*/
        /*std::hypot(resolution.second, 1e2 * shifts.at(staveId).second)};*/
        resolution.first, resolution.second};
    Acts::SquareMatrix2 cov = stdDev.cwiseProduct(stdDev).asDiagonal();
    return {cov, digLocal};
  }
};
