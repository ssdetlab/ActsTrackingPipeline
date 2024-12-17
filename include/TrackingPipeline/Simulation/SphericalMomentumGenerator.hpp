#pragma once

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"

/// @brief Spherical momentum generator
struct SphericalMomentumGenerator : public IMomentumGenerator {
    std::pair<Acts::ActsScalar, Acts::ActsScalar> pRange;
    std::pair<Acts::ActsScalar, Acts::ActsScalar> thetaRange;
    std::pair<Acts::ActsScalar, Acts::ActsScalar> phiRange;

    Acts::Vector3 gen(RandomEngine& rng) const override {
        std::uniform_real_distribution<Acts::ActsScalar> uniform(0, 1);
        
        Acts::ActsScalar pMag = 
            pRange.first + (pRange.second - pRange.first) * uniform(rng);

        Acts::ActsScalar phi = 
            phiRange.first + (phiRange.second - phiRange.first) * uniform(rng);

        Acts::ActsScalar theta = 
            std::acos(
                std::cos(thetaRange.first) - 
                    (std::cos(thetaRange.first) - std::cos(thetaRange.second)) * 
                        uniform(rng));

        return pMag * Acts::Vector3(
            std::sin(theta) * std::cos(phi),
            std::sin(theta) * std::sin(phi),
            std::cos(theta));
    }
};
