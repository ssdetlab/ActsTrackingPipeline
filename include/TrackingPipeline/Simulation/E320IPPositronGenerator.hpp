#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <memory>
#include <vector>

#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Io/ITrackParamsReader.hpp"
#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

namespace E320Sim {

class E320IPPositronGenerator : public IVertexGenerator,
                                public IMomentumGenerator {
 public:
  struct Config {
    std::shared_ptr<ITrackParamsReader> trackParamsReader;
  };

  struct State {
    Acts::Vector3 vertex;
    Acts::Vector3 momentum;

    std::size_t idx;

    std::pair<bool, bool> genState;
  };

  E320IPPositronGenerator(const Config& cfg)
      : m_cfg(cfg),
        m_params(m_cfg.trackParamsReader->read()),
        m_state(std::make_unique<State>()) {
    m_state->genState = {false, false};
    m_state->idx = 0;
  };

  Acts::Vector3 genVertex(RandomEngine& rng) const override {
    // std::cout << "\n\n\n\n";
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

  std::vector<Acts::CurvilinearTrackParameters> m_params;

  void internalUpdate(RandomEngine& /*rng*/) const {
    m_state->vertex = m_params.at(m_state->idx).position();
    m_state->momentum = m_params.at(m_state->idx).momentum();
    m_state->genState = {true, true};
    m_state->idx++;
  }

  E320Geometry::GeometryOptions m_gOpt;
};

}  // namespace E320Sim
