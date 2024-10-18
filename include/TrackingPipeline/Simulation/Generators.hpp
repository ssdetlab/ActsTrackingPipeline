#pragma once

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/LUXEGeometryConstraints.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

/// @brief Interface for generating vertex positions
struct IVertexGenerator {
    virtual Acts::Vector3 gen(RandomEngine& rng) const = 0;
};

/// @brief Uniform vertex generator
struct UniformVertexGenerator : public IVertexGenerator {
    Acts::Vector3 mins{0., 0., 0.};
    Acts::Vector3 maxs{0., 0., 0.};

    Acts::Vector3 gen(RandomEngine& rng) const override {
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

    Acts::Vector3 gen(RandomEngine& /*rng*/) const override {
        return vertex;
    }
};

/// @brief Interface for generating momentum vectors
struct IMomentumGenerator {
    virtual Acts::Vector3 gen(RandomEngine& rng) const = 0;
};

/// @brief Gaussian momentum generator
struct GaussianMomentumGenerator : public IMomentumGenerator {
    std::pair<Acts::ActsScalar, Acts::ActsScalar> pMagRange;
    std::pair<Acts::ActsScalar, Acts::ActsScalar> thetaRange;
    std::pair<Acts::ActsScalar, Acts::ActsScalar> phiRange;

    Acts::Vector3 gen(RandomEngine& rng) const override {
        std::uniform_real_distribution<Acts::ActsScalar> uniform(0, 1);

        Acts::ActsScalar pMag = pMagRange.first + (pMagRange.second - pMagRange.first) * uniform(rng);

        Acts::ActsScalar theta = std::acos(
            1 - (std::cos(thetaRange.first) - std::cos(thetaRange.second)) * uniform(rng));

        Acts::ActsScalar phi = phiRange.first + (phiRange.second - phiRange.first) * uniform(rng);

        return pMag * Acts::Vector3(
            std::sin(theta) * std::cos(phi),
            std::cos(theta),
            -std::sin(theta) * std::sin(phi));
    }
};

/// @brief Uniform momentum generator
struct RangedUniformMomentumGenerator : public IMomentumGenerator {
    std::vector<std::pair<Acts::ActsScalar, Acts::ActsScalar>> Pranges;

    Acts::Vector3 gen(RandomEngine& rng) const override {
        std::uniform_int_distribution<int> range_select(0, Pranges.size() - 1);
        int range = range_select(rng);

        Acts::ActsScalar Pmin = Pranges.at(range).first;
        Acts::ActsScalar Pmax = Pranges.at(range).second;

        std::uniform_real_distribution<Acts::ActsScalar> uniform(Pmin, Pmax);
        Acts::ActsScalar p = uniform(rng);

        return p * Acts::Vector3(0, 1, 0);
    }
};

/// @brief Interface for generating random noise hits
struct INoiseGenerator {
    virtual std::vector<SimpleSourceLink> gen(
        RandomEngine& rng, const Acts::Surface* surface) const = 0;
};

/// @brief Uniform noise generator
struct UniformNoiseGenerator : public INoiseGenerator {
    Acts::ActsScalar numberOfHits = 10;

    std::vector<SimpleSourceLink> gen(
        RandomEngine& rng, const Acts::Surface* surface) const override {
            if (surface->type() != Acts::Surface::SurfaceType::Plane) {
                throw std::runtime_error("Only plane surfaces are supported");
            }
            if (surface->geometryId().sensitive() == 0) {
                return {};
            }
            std::uniform_real_distribution<Acts::ActsScalar> uniform(0, 1);
    
            auto surfaceBounds = surface->bounds().values();
            Acts::ActsScalar area = 
                (surfaceBounds.at(3) - surfaceBounds.at(1)) * 
                (surfaceBounds.at(2) - surfaceBounds.at(0));
            
            std::vector<SimpleSourceLink> noiseHits;    
            while (noiseHits.size() < numberOfHits) {
                Acts::ActsScalar x = 
                    surfaceBounds.at(0) + (surfaceBounds.at(2) - surfaceBounds.at(0)) * uniform(rng);
                Acts::ActsScalar y = 
                    surfaceBounds.at(1) + (surfaceBounds.at(3) - surfaceBounds.at(1)) * uniform(rng);
                Acts::Vector2 loc{x, y};

                Acts::Vector2 stddev(5  * Acts::UnitConstants::um,
                    5  * Acts::UnitConstants::um);
                Acts::SquareMatrix2 cov = stddev.cwiseProduct(stddev).asDiagonal();
                SimpleSourceLink simpleSL(loc, cov, surface->geometryId(), -1, -1);
                noiseHits.push_back(simpleSL);
            }

            return noiseHits;
    }
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

} // namespace LUXESimParticle
