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
        
        Acts::ActsScalar meanP = 2_MeV;
        Acts::ActsScalar sigmaP = 1.8_MeV;

        Acts::ActsScalar meanM = -2_MeV;
        Acts::ActsScalar sigmaM = 1.8_MeV;

        Acts::ActsScalar shapeZ = 3_GeV;
        Acts::ActsScalar scaleZ = 1.2_GeV;

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

    /// @brief Uniform momentum generator
    struct RangedUniformMomentumGenerator : public IMomentumGenerator {
        std::vector<std::pair<Acts::ActsScalar, Acts::ActsScalar>> Eranges {
            std::make_pair(1_GeV, 3_GeV),
            std::make_pair(3_GeV, 6_GeV),
            std::make_pair(6_GeV, 9_GeV),
            std::make_pair(9_GeV, 12_GeV)};
    
        Acts::ActsScalar m = 0.511 * Acts::UnitConstants::MeV;
    
        Acts::Vector3 gen(RandomEngine rng) const override {
            std::uniform_int_distribution<int> range_select(0, Eranges.size() - 1);
            int range = range_select(rng);
    
            Acts::ActsScalar Emin = Eranges.at(range).first;
            Acts::ActsScalar Emax = Eranges.at(range).second;
    
            std::uniform_real_distribution<Acts::ActsScalar> uniform(Emin, Emax);
            Acts::ActsScalar E = uniform(rng);
            Acts::ActsScalar p = std::sqrt(E * E - m * m);
    
            return p * Acts::Vector3(0, 1, 0);
        }
    };
}
