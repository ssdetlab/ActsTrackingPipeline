#pragma once

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

#include "Acts/Definitions/Algebra.hpp"

/// @brief Class that digitizes hits based on the provided 
/// resolution assuming Gaussian smearing
struct SimpleDigitizer : public IDigitizer {
    std::pair<Acts::ActsScalar, Acts::ActsScalar> resolution;

    std::pair<Acts::SquareMatrix2, Acts::Vector2> gen(
        RandomEngine& rng,
        Acts::GeometryIdentifier /*geoId*/,
        Acts::Vector2 pos) const override {
            std::normal_distribution<double> normalDist(0., 1.);
    
            Acts::Vector2 stdDev = {resolution.first, resolution.second};
            Acts::Vector2 digLocal = pos + stdDev.cwiseProduct(
                Acts::Vector2(normalDist(rng), normalDist(rng)));
    
            Acts::SquareMatrix2 cov = stdDev.cwiseProduct(stdDev).asDiagonal();
    
            return {cov, digLocal};
    }
};
