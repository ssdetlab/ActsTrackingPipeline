#pragma once

#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Io/ITrackParamsReader.hpp"
#include "TrackingPipeline/Simulation/PowerLawVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/detail/NormalKDE.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <memory>
#include <vector>

namespace E320Sim {

class E320CptBkgGenerator : public IVertexGenerator, public IMomentumGenerator {
    public: 
        struct Config {
            std::shared_ptr<ITrackParamsReader> trackParamsReader;
            std::size_t nIterations;
            Acts::ActsScalar sensitivity;
            Acts::ActsScalar yPower;
            Acts::ActsScalar yShift;
            std::vector<Acts::ActsScalar> zProbs;
        };

        struct State {
            Acts::Vector3 vertex;
            Acts::Vector3 momentum;
            
            std::pair<bool, bool> genState;
        };

        E320CptBkgGenerator(const Config& cfg) 
            : m_cfg(cfg),
            m_state(std::make_unique<State>()),
            m_vertexGen(std::make_unique<E320PowerLawVertexGenerator>()) {
                m_state->genState = {false, false};
    
                auto trackParams = m_cfg.trackParamsReader->read();
    
                std::vector<Acts::Vector4> sample;
                sample.reserve(trackParams.size());
                for (const auto& param : trackParams) {
                    Acts::ActsScalar zSub = 0;
                    if (std::abs(param.fourPosition().z() - 
                        m_gOpt.staveZ.at(0)) < 10) {
                            zSub = 0;
                    }
                    else if (std::abs(param.fourPosition().z() - 
                        m_gOpt.staveZ.at(1)) < 10) {
                            zSub = 1;
                    }
                    else if (std::abs(param.fourPosition().z() - 
                        m_gOpt.staveZ.at(2)) < 10) {
                            zSub = 2;
                    }
                    else if (std::abs(param.fourPosition().z() - 
                        m_gOpt.staveZ.at(3)) < 10) {
                            zSub = 3;
                    }
    
                    sample.emplace_back(
                        Acts::Vector4{
                            zSub,
                            param.phi(), 
                            param.theta(),
                            param.absoluteMomentum()});
                }
                m_zPhiThetaEKDE = std::make_unique<NormalKDE<4>>(
                    std::move(sample),
                    m_cfg.nIterations,
                    m_cfg.sensitivity);
    
                m_vertexGen->yBoundLow = m_gOpt.chipY.at(0) - m_gOpt.chipSizeY/2;
                m_vertexGen->yBoundHigh = m_gOpt.chipY.at(8) + m_gOpt.chipSizeY/2;
        
                m_vertexGen->xBoundLow = m_gOpt.chipX - m_gOpt.chipSizeX/2;
                m_vertexGen->xBoundHigh = m_gOpt.chipX + m_gOpt.chipSizeX/2;

                m_vertexGen->yPower = m_cfg.yPower;
                m_vertexGen->yShift = m_cfg.yShift;
                m_vertexGen->zProbs = m_cfg.zProbs;
                m_vertexGen->zPositions = {
                    m_gOpt.staveZ.at(0),
                    m_gOpt.staveZ.at(1),
                    m_gOpt.staveZ.at(2),
                    m_gOpt.staveZ.at(3)};
        };

        Acts::Vector3 genVertex(RandomEngine& rng) const override {
            if (!m_state->genState.first) {
                internalUpdate(rng);
            }
            m_state->genState.first = false;
            return m_state->vertex;
        };

        Acts::Vector3 genMomentum(RandomEngine& rng) const override {
            if (!m_state->genState.second) {
                internalUpdate(rng);
            }
            m_state->genState.second = false;
            return m_state->momentum;
        };

    private:
        Config m_cfg;

        std::unique_ptr<State> m_state;

        std::unique_ptr<NormalKDE<4>> m_zPhiThetaEKDE;

        std::unique_ptr<E320PowerLawVertexGenerator> m_vertexGen;

        void internalUpdate(RandomEngine &rng) const {
            Acts::Vector3 vertex = m_vertexGen->genVertex(rng);
            Acts::Vector4 zPhiThetaE = m_zPhiThetaEKDE->sample(rng);

            int zIdx = static_cast<int>(std::round(zPhiThetaE(0)));

            Acts::ActsScalar phi = zPhiThetaE(1);
            Acts::ActsScalar theta = zPhiThetaE(2);
            Acts::ActsScalar E = zPhiThetaE(3);

            Acts::Vector3 dir{
                std::sin(theta) * std::cos(phi),
                std::sin(theta) * std::sin(phi),
                std::cos(theta)};
            dir = m_gOpt.actsToWorld.rotation().inverse() * dir;

            Acts::ActsScalar z = m_gOpt.staveZ.at(zIdx);

            Acts::Vector3 pos = Acts::Vector3{
                vertex(0),
                z,
                vertex(2)};

            m_state->vertex = pos;
            m_state->momentum = E * dir;
            m_state->genState = {true, true};
        }

        E320Geometry::GeometryOptions m_gOpt;
};

} // namespace E320Sim