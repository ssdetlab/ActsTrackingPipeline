#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <random>
#include <stdexcept>
#include <vector>

#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

namespace E320Sim {

/// @brief Vertex sampler constructed to sample the inital
/// positions of the NCS-induced backgorund in the tracking
/// detector of the E320 experiment
class E320CptBkgVertexGenerator : public IVertexGenerator {
 public:
  struct Config {
    /// y-powerlaw parameters
    double yBoundLow;
    double yBoundHigh;
    double yScale;
    double yShift;
    double yPedestal;
    double yNorm;
    double xBoundLow;
    double xBoundHigh;
    std::vector<double> zProbs;
    std::vector<double> zPositions;
  };

  E320CptBkgVertexGenerator(const Config& cfg)
      : m_yBoundLow(cfg.yBoundLow),
        m_yBoundHigh(cfg.yBoundHigh),
        m_yScale(cfg.yScale),
        m_yShift(cfg.yShift),
        m_yPedestal(cfg.yPedestal),
        m_yNorm(cfg.yNorm),
        m_xBoundLow(cfg.xBoundLow),
        m_xBoundHigh(cfg.xBoundHigh),
        m_zProbs(cfg.zProbs),
        m_zPositions(cfg.zPositions),
        m_is({m_yBoundLow, m_yBoundHigh}),
        m_envelopeA((powerLaw(m_yBoundHigh) - powerLaw(m_yBoundLow)) /
                    (m_yBoundHigh - m_yBoundLow)),
        m_envelopeB(powerLaw(m_yBoundLow) - m_envelopeA * m_yBoundLow + 1e-4),
        m_ws({powerLaw(m_yBoundLow), powerLaw(m_yBoundHigh)}) {}

  Acts::Vector3 genVertex(RandomEngine& rng) const override {
    if (m_zProbs.empty()) {
      throw std::runtime_error("Vertex generator not initialized");
    }

    std::discrete_distribution<> discrete(m_zProbs.begin(), m_zProbs.end());
    std::uniform_real_distribution<> uniform(0, 1);

    std::piecewise_linear_distribution<> linear(m_is.begin(), m_is.end(),
                                                m_ws.begin());

    // Generate unfiorm x
    double x = m_xBoundLow + (m_xBoundHigh - m_xBoundLow) * uniform(rng);

    // Generate power-law y
    bool sampled = false;
    double y;
    while (!sampled) {
      // Sample envelope
      double envY = linear(rng);

      // Sample ratio
      double ratio = uniform(rng);

      // Check ratio
      if (ratio < powerLaw(envY) / (m_envelopeA * envY + m_envelopeB)) {
        y = envY;
        sampled = true;
      }
    }

    // Generate uniform z
    double z = m_zPositions.at(discrete(rng));

    Acts::Vector3 glob(x, z, -y);
    return glob;
  }

 private:
  double m_yBoundLow;
  double m_yBoundHigh;

  double m_yScale;
  double m_yShift;
  double m_yPedestal;
  double m_yNorm;

  double m_xBoundLow;
  double m_xBoundHigh;

  std::vector<double> m_zProbs;
  std::vector<double> m_zPositions;

  std::vector<double> m_is;
  std::vector<double> m_ws;
  double m_envelopeA;
  double m_envelopeB;

  double powerLaw(double x) const {
    return (m_yScale * std::pow(x + m_yShift, -0.01) + m_yPedestal) / m_yNorm;
  }
};

inline E320CptBkgVertexGenerator::Config cptBkgVertexGenConfig() {
  E320CptBkgVertexGenerator::Config cfg;
  E320Geometry::GeometryOptions gOpt;
  // y-power law
  cfg.yBoundLow = gOpt.chipY.at(0) - gOpt.chipSizeY / 2;
  cfg.yBoundHigh = gOpt.chipY.at(8) + gOpt.chipSizeY / 2;
  cfg.yScale = 5.31456e3;
  cfg.yShift = -6.98928e1;
  cfg.yPedestal = -4.99682e3;
  cfg.yNorm = 20011.3;
  // x-uniform
  cfg.xBoundLow = gOpt.chipX - gOpt.chipSizeX;
  cfg.xBoundHigh = gOpt.chipX + gOpt.chipSizeX;
  // z-discrete
  cfg.zProbs = {0.25, 0.24, 0.23, 0.28};
  cfg.zPositions = {gOpt.staveZ.at(0), gOpt.staveZ.at(1), gOpt.staveZ.at(2),
                    gOpt.staveZ.at(3)};
  return cfg;
}

}  // namespace E320Sim
