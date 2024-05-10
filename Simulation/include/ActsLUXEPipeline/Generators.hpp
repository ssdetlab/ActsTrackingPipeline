#pragma once

#include "ActsLUXEPipeline/RandomNumbers.hpp"
#include"ActsLUXEPipeline/LUXEGeometryConstraints.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"

/// @brief Interface for generating vertex positions
struct IVertexGenerator {
    virtual Acts::Vector3 gen(RandomEngine rng) const = 0;
};

/// @brief Uniform vertex generator
struct UniformVertexGenerator : public IVertexGenerator {
    Acts::Vector3 mins{0., 0., 0.};
    Acts::Vector3 maxs{0., 0., 0.};

    Acts::Vector3 gen(RandomEngine rng) const override {
        std::uniform_real_distribution<Acts::ActsScalar> uniform;
        Acts::Vector3 vertex{
            uniform(rng),
            uniform(rng),
            uniform(rng)};
        return mins + vertex.cwiseProduct(maxs - mins);
    }
};

/// @brief Stationary vertex generator
struct StationaryVertexGenerator : public IVertexGenerator {
    Acts::Vector3 vertex{0., 0., 0.};

    Acts::Vector3 gen(RandomEngine /*rng*/) const override {
        return vertex;
    }
};

/// @brief Interface for generating momentum vectors
struct IMomentumGenerator {
    virtual Acts::Vector3 gen(RandomEngine rng) const = 0;
};

/// @brief LUXE specific generators
namespace LUXESimParticle {
    using namespace Acts::UnitLiterals;

    struct LUXEMomentumGenerator : public IMomentumGenerator {
        LUXEGeometry::GeometryOptions gOpt;
        
        Acts::ActsScalar meanP = 2 * 1_MeV;
        Acts::ActsScalar sigmaP = 1.8 * 1_MeV;

        Acts::ActsScalar meanM = -2 * 1_MeV;
        Acts::ActsScalar sigmaM = 1.8 * 1_MeV;

        Acts::ActsScalar shapeZ = 3;
        Acts::ActsScalar scaleZ = 1.2;


        Acts::Vector3 gen(RandomEngine rng) const override {
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

}
