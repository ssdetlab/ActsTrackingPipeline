#pragma once

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/detail/NormalRandomVariable.hpp"

namespace E320Sim {

    struct E320MomentumGenerator : public IMomentumGenerator {
        NormalRandomVariable normal;

        Acts::ActsScalar gmDistMean;
        Acts::ActsScalar gmDistSigma;
        Acts::ActsScalar gmDistScale;

        Acts::ActsScalar gvDistMean;
        Acts::ActsScalar gvDistSigma;
        Acts::ActsScalar gvDistScale;

        E320MomentumGenerator(
            Acts::Vector2 mean, 
            Acts::SquareMatrix2 cov,
            Acts::ActsScalar gmMean,
            Acts::ActsScalar gmSigma,
            Acts::ActsScalar gmScale,
            Acts::ActsScalar gvMean,
            Acts::ActsScalar gvSigma,
            Acts::ActsScalar gvScale)
                : normal(mean, cov),
                gmDistMean(gmMean),
                gmDistSigma(gmSigma),
                gmDistScale(gmScale),
                gvDistMean(gvMean),
                gvDistSigma(gvSigma),
                gvDistScale(gvScale){};

        Acts::ActsScalar normalPDF(
            Acts::ActsScalar x, 
            Acts::ActsScalar mu, 
            Acts::ActsScalar sigma, 
            Acts::ActsScalar scale) const {
                return scale / (sigma * std::sqrt(2 * M_PI)) * 
                    std::exp(-0.5 * (x - mu) * (x - mu) / (sigma * sigma));
        }

        Acts::Vector3 gen(RandomEngine& rng) const override {
            Acts::Vector2 angles = normal.gen(rng);

            Acts::ActsScalar theta = angles.x();
            Acts::ActsScalar phi = angles.y();

            double gammaMean = normalPDF(
                phi, 
                gmDistMean, 
                gmDistSigma, 
                gmDistScale);
            double gammaVar = normalPDF(
                phi, 
                gvDistMean, 
                gvDistSigma, 
                gvDistScale);
    
            double beta = gammaVar / gammaMean;
            double alpha = gammaMean / beta;

            std::gamma_distribution<> distE(alpha, beta);
            Acts::ActsScalar E = distE(rng);

            Acts::Vector3 dir{
                std::sin(theta) * std::cos(phi),
                std::sin(theta) * std::sin(phi),
                std::cos(theta)};

            return E * dir;
        }
    };

} // namespace E320Sim
