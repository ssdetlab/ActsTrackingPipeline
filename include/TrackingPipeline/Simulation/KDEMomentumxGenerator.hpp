#pragma once

#include "TrackingPipeline/Io/ITrackParamsReader.hpp"

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/detail/NormalKDE.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <memory>
#include <vector>

/// @brief Uniform momentum generator
class KDEMomentumGenerator : public IMomentumGenerator {
    public:
        struct Config {
            std::shared_ptr<ITrackParamsReader> trackParamsReader;
            std::size_t nIterations;
            Acts::ActsScalar sensitivity;
            Acts::Transform3 transform;
        };
    
        KDEMomentumGenerator(const Config& cfg) : m_cfg(cfg) {
            auto trackParams = m_cfg.trackParamsReader->read();
    
            std::vector<Acts::Vector3> sample;
            sample.reserve(trackParams.size());
            for (const auto& param : trackParams) {
                sample.emplace_back(
                    Acts::Vector3{
                        param.phi(), 
                        param.theta(), 
                        param.absoluteMomentum()});
            }
            m_kde = std::make_unique<NormalKDE<3>>(
                std::move(sample),
                m_cfg.nIterations,
                m_cfg.sensitivity);
        }

        Acts::Vector3 genMomentum(RandomEngine& rng) const override {
            Acts::Vector3 phiThetaE = m_kde->sample(rng);
            Acts::ActsScalar phi = phiThetaE.x();
            Acts::ActsScalar theta = phiThetaE.y();
            Acts::ActsScalar E = phiThetaE.z();

            Acts::Vector3 dir{
                std::sin(theta) * std::cos(phi),
                std::sin(theta) * std::sin(phi),
                std::cos(theta)};

            dir = m_cfg.transform * dir;

            return E * dir;
        }
    
    private:
        Config m_cfg;

        std::unique_ptr<NormalKDE<3>> m_kde;
};
