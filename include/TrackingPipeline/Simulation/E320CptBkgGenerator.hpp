#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <memory>
#include <vector>

#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Io/ITrackParamsReader.hpp"
#include "TrackingPipeline/Simulation/E320CptBkgVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/detail/NormalKDE.hpp"

namespace E320Sim {

class E320CptBkgGenerator : public IVertexGenerator, public IMomentumGenerator {
 public:
  struct Config {
    std::shared_ptr<ITrackParamsReader> trackParamsReader;
    std::size_t nIterations;
    double sensitivity;
    E320CptBkgVertexGenerator::Config vertexGenCfg;
  };

  struct State {
    Acts::Vector3 vertex;
    Acts::Vector3 momentum;

    std::pair<bool, bool> genState;
  };

  E320CptBkgGenerator(const Config& cfg)
      : m_cfg(cfg),
        m_state(std::make_unique<State>()),
        m_vertexGen(
            std::make_unique<E320CptBkgVertexGenerator>(m_cfg.vertexGenCfg)) {
    m_state->genState = {false, false};

    auto trackParams = m_cfg.trackParamsReader->read();

    std::vector<Acts::Vector4> sample;
    sample.reserve(trackParams.size());
    for (const auto& param : trackParams) {
      double zSub = 0;
      if (std::abs(param.fourPosition().z() - m_gOpt.staveZ.at(0)) < 10) {
        zSub = 0;
      } else if (std::abs(param.fourPosition().z() - m_gOpt.staveZ.at(1)) <
                 10) {
        zSub = 1;
      } else if (std::abs(param.fourPosition().z() - m_gOpt.staveZ.at(2)) <
                 10) {
        zSub = 2;
      } else if (std::abs(param.fourPosition().z() - m_gOpt.staveZ.at(3)) <
                 10) {
        zSub = 3;
      }

      sample.emplace_back(Acts::Vector4{zSub, param.phi(), param.theta(),
                                        param.absoluteMomentum()});
    }

    m_zPhiThetaEKDE = std::make_unique<NormalKDE<4>>(
        std::move(sample), m_cfg.nIterations, m_cfg.sensitivity);
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

  std::unique_ptr<E320CptBkgVertexGenerator> m_vertexGen;

  void internalUpdate(RandomEngine& rng) const {
    Acts::Vector3 vertex = m_vertexGen->genVertex(rng);
    Acts::Vector4 zPhiThetaE = m_zPhiThetaEKDE->sample(rng);

    int zIdx = static_cast<int>(std::round(zPhiThetaE(0)));

    double phi = zPhiThetaE(1);
    double theta = zPhiThetaE(2);
    double E = zPhiThetaE(3);

    Acts::Vector3 dir{std::sin(theta) * std::cos(phi),
                      std::sin(theta) * std::sin(phi), std::cos(theta)};
    dir = m_gOpt.actsToWorld.rotation().inverse() * dir;

    double z = m_gOpt.staveZ.at(zIdx);

    Acts::Vector3 pos = Acts::Vector3{vertex(0), z, vertex(2)};

    m_state->vertex = pos;
    m_state->momentum = E * dir;
    m_state->genState = {true, true};
  }

  E320Geometry::GeometryOptions m_gOpt;
};

}  // namespace E320Sim
