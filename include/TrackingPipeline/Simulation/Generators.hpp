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
        const Acts::GeometryContext& geoContext,
        RandomEngine& rng, 
        const std::vector<const Acts::Surface*> surfaces) const = 0;
};

/// @brief Uniform noise generator
struct UniformNoiseGenerator : public INoiseGenerator {
    Acts::ActsScalar numberOfHits = 10;

    std::vector<SimpleSourceLink> gen(
        const Acts::GeometryContext& /*geoContext*/,
        RandomEngine& rng, 
        const std::vector<const Acts::Surface*> surfaces) const override {
            std::uniform_real_distribution<Acts::ActsScalar> uniform(0, 1);

            std::vector<SimpleSourceLink> noiseHits;    
            for (auto surface : surfaces) {
                if (surface->type() != Acts::Surface::SurfaceType::Plane) {
                    throw std::runtime_error("Only plane surfaces are supported");
                }
                if (surface->geometryId().sensitive() == 0) {
                    return {};
                }
                auto surfaceBounds = surface->bounds().values();
                Acts::ActsScalar area = 
                    (surfaceBounds.at(3) - surfaceBounds.at(1)) * 
                    (surfaceBounds.at(2) - surfaceBounds.at(0));
                
                for (int i = 0; i < numberOfHits; i++) {
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
            }
            return noiseHits;
    }
};

/// @brief E320 specific generators
namespace E320Sim {
    using namespace Acts::UnitLiterals;

    struct E320ComptonBackgroundGenerator : public INoiseGenerator {
        /// Number of background hits
        /// in an event
        int numberOfHits = 10;

        /// Geometry constraints
        E320Geometry::GeometryOptions gOpt;

        /// Boundary tolerance
        Acts::BoundaryTolerance boundaryTolerance = 
            Acts::BoundaryTolerance::None();

        /// y power-law parameters
        double yPower = -2.66;
        double yShift = -158.0;
        double yMin = gOpt.chipY.at(0) - gOpt.chipSizeY / 2;
        double yMax = gOpt.chipY.at(8) + gOpt.chipSizeY / 2;
    
        /// x uniform parameters
        double xMin = gOpt.chipX - gOpt.chipSizeX / 2;
        double xMax = gOpt.chipX + gOpt.chipSizeX / 2;
    
        /// Size 2D discrete distribution
        /// with mapping (SizeX, SizeY) -> SizeX * 12 + SizeY
        /// as the 12 is going to be the maximum number of
        /// pixels in x/y direction
        std::array<double, 144> sizeProbs;
    
        /// Table with sparse probabilities
        /// for the size distribution
        std::map<std::pair<int, int>, double> sizeProbTable = {
            {{2, 2}, 0.121},
            {{2, 3}, 0.141},
            {{2, 4}, 0.000417},
            {{3, 2}, 0.082},
            {{3, 3}, 0.457},
            {{3, 4}, 0.0777},
            {{3, 5}, 0.00177},
            {{3, 6}, 0.000209},
            {{4, 3}, 0.0248},
            {{4, 4}, 0.0455},
            {{4, 5}, 0.0148},
            {{4, 6}, 0.00115},
            {{4, 7}, 0.000522},
            {{4, 8}, 0.000104},
            {{4, 9}, 0.000209},
            {{4, 10}, 0.000104},
            {{5, 3}, 0.000313},
            {{5, 4}, 0.00678},
            {{5, 5}, 0.0135},
            {{5, 6}, 0.00261},
            {{5, 7}, 0.000939},
            {{5, 8}, 0.000209},
            {{5, 9}, 0.000209},
            {{6, 3}, 0.000313},
            {{6, 4}, 0.000417},
            {{6, 5}, 0.0023},
            {{6, 6}, 0.00125},
            {{6, 7}, 0.000522},
            {{6, 8}, 0.000209},
            {{6, 9}, 0.000104},
            {{7, 4}, 0.000104},
            {{7, 5}, 0.000104},
            {{7, 6}, 0.000313},
            {{7, 7}, 0.000209},
            {{7, 8}, 0.000104},
            {{7, 9}, 0.000104},
            {{8, 4}, 0.000104},
            {{8, 5}, 0.000209},
            {{8, 6}, 0.000104},
            {{9, 6}, 0.000209},
            {{9, 8}, 0.000104},
            {{10, 9}, 0.000209},
            {{10, 10}, 0.000104}
        };

        E320ComptonBackgroundGenerator() {
            // Fill the size distribution
            for (int i = 1; i <= 12; i++) {
                for (int j = 1; j <= 12; j++) {
                    // Some sizes are not present in the table
                    if (sizeProbTable.find({i, j}) == sizeProbTable.end()) {
                        sizeProbs.at(
                            (i - 1) * 12 + (j - 1)) = 0;
                    } else {
                        sizeProbs.at(
                            (i - 1) * 12 + (j - 1)) = sizeProbTable.at({i, j});
                    }
                }
            }
        }

        std::vector<SimpleSourceLink> gen(
            const Acts::GeometryContext& geoContext,
            RandomEngine& rng, 
            const std::vector<const Acts::Surface*> surfaces) const override {
                for (auto surface : surfaces) {
                    if (surface->type() != Acts::Surface::SurfaceType::Plane) {
                        throw std::runtime_error("Only plane surfaces are supported");
                    }
                }

                // For error estimation
                double pixSizeX = 27_um;
                double pixSizeY = 29_um;

                // Generate uniform x
                std::uniform_real_distribution<double> uniform(0, 1);

                std::vector<SimpleSourceLink> noiseHits;    
                for (int i = 0; i < numberOfHits; i++) {
                    // Generate x uniform
                    double x = xMin + (xMax - xMin) * uniform(rng);
    
                    // Generate y power law
                    // via the inverse transform sampling
                    double gamma = uniform(rng);

                    double term1 = gamma * (
                        std::pow(yMax - yShift, yPower + 1) - 
                        std::pow(yMin - yShift, yPower + 1));
            
                    double term2 = std::pow(yMin - yShift, yPower + 1);
            
                    double y = yShift + std::pow(term1 + term2, 1 / (yPower + 1));
            
                    // Generate uniform z
                    std::discrete_distribution<int> zDist({
                        0.25, 0.25, 0.25, 0.25});
            
                    double z = 0;
                    switch (zDist(rng))
                    {
                        case 0:
                            z = 16567.0_mm;
                            break;
                        case 1:
                            z = 16667.0_mm;
                            break;
                        case 2:
                            z = 16767.0_mm;
                            break;
                        case 3:
                            z = 16867.0_mm;
                            break;
                        default:
                            break;
                    }
            
                    // Generate correlated sizes
                    std::discrete_distribution<int> sizeDist(
                        sizeProbs.begin(), sizeProbs.end());
            
                    int idx = sizeDist(rng);
            
                    int sizeX = idx / 12;
                    int sizeY = idx % 12;

                    Acts::Vector3 glob(x, y, z);
                    glob = gOpt.actsToWorld.rotation().inverse() * glob;

                    Acts::Vector2 stddev(sizeX * pixSizeX, sizeY * pixSizeY);
                    Acts::SquareMatrix2 cov = stddev.cwiseProduct(stddev).asDiagonal();

                    for (auto surface : surfaces) {
                        Acts::Vector3 center = surface->center(geoContext);

                        if (surface->isOnSurface(
                                geoContext, 
                                glob, 
                                Acts::Vector3(0, 1, 0), 
                                boundaryTolerance, 
                                1e-1)) {
                                    Acts::Vector2 loc = surface->globalToLocal(
                                        geoContext, glob, Acts::Vector3(0, 1, 0), 1e-1).value();

                                    SimpleSourceLink simpleSL(
                                        loc, 
                                        cov, 
                                        surface->geometryId(), 
                                        -1, 
                                        -1);
                                    noiseHits.push_back(simpleSL);
                                    break;
                        }
                    }
                }

                return noiseHits;
        }
    };


} // namespace E320Sim

/// @brief LUXE specific generators
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
