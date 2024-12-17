#pragma once

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Geometry/LUXEGeometryConstraints.hpp"

#include "Acts/Definitions/Units.hpp"

namespace LUXESim {

    using namespace Acts::UnitLiterals;

    struct LUXEMomentumGenerator : public IMomentumGenerator {
        LUXEGeometry::GeometryOptions gOpt;
        
        Acts::ActsScalar meanP = 2_MeV;
        Acts::ActsScalar sigmaP = 1.8_MeV;

        Acts::ActsScalar meanM = -2_MeV;
        Acts::ActsScalar sigmaM = 1.8_MeV;

        Acts::ActsScalar shapeZ = 3_GeV;
        Acts::ActsScalar scaleZ = 1.2_GeV;

        Acts::Vector3 gen(RandomEngine& rng) const override {
            std::normal_distribution<Acts::ActsScalar> pDisP(meanP, sigmaP);
            std::normal_distribution<Acts::ActsScalar> pDisM(meanM, sigmaM);
            std::gamma_distribution<Acts::ActsScalar> pzDis(shapeZ, scaleZ);

            Acts::ActsScalar px = (pDisP(rng) + pDisM(rng))/2;
            Acts::ActsScalar py = (pDisP(rng) + pDisM(rng))/2;
            Acts::ActsScalar pz = pzDis(rng) + 1;

            Acts::Vector3 momentum{px, py, pz};

            auto worldMomentum = gOpt.actsToWorld.rotation().inverse() * momentum;
            return worldMomentum;
        }
    };

} // namespace LUXESim
