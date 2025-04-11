#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <random>
#include <stdexcept>
#include <vector>

#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

namespace E320Sim {

class E320CptBkgVertexGenerator : public IVertexGenerator {
 public:
  struct Config {
    double yBoundLow;
    double yBoundHigh;

    double xBoundLow;
    double xBoundHigh;

    double yScale;
    double yShift;
    double yPedestal;
    double yPower;
    double yNorm;

    std::vector<double> zProbs;
    std::vector<double> zPositions;
  };

  E320CptBkgVertexGenerator(const Config& cfg)
      : m_cfg(cfg),
        m_yBoundLow(m_cfg.yBoundLow),
        m_yBoundHigh(m_cfg.yBoundHigh),
        m_xBoundLow(m_cfg.xBoundLow),
        m_xBoundHigh(m_cfg.xBoundHigh),
        m_yScale(m_cfg.yScale),
        m_yShift(m_cfg.yShift),
        m_yPedestal(m_cfg.yPedestal),
        m_yPower(m_cfg.yPower),
        m_yNorm(m_cfg.yNorm),
        m_zProbs(m_cfg.zProbs),
        m_zPositions(m_cfg.zPositions) {
    m_is = std::vector<double>{m_yBoundLow, m_yBoundHigh};
    m_ws = std::vector<double>{powerLaw(m_yBoundLow), powerLaw(m_yBoundHigh)};
    m_linearA = (powerLaw(m_yBoundHigh) - powerLaw(m_yBoundLow)) /
                (m_yBoundHigh - m_yBoundLow);
    m_linearB = powerLaw(m_yBoundLow) - m_linearA * m_yBoundLow + 1e-4;
  };

  Acts::Vector3 genVertex(RandomEngine& rng) const override {
    if (m_zProbs.size() == 0) {
      throw std::runtime_error("Vertex generator not initialized");
    }

    std::discrete_distribution<> discrete(m_zProbs.begin(), m_zProbs.end());
    std::uniform_real_distribution<> uniform(0, 1);

    std::piecewise_linear_distribution<> linear(m_is.begin(), m_is.end(),
                                                m_ws.begin());

    // Generate unfiorm x
    double x = m_xBoundLow + (m_xBoundHigh - m_xBoundLow) * uniform(rng);

    // Generate power-law y
    double y;
    bool sampled = false;
    while (!sampled) {
      // Sample envelope
      double envY = linear(rng);

      // Sample ratio
      double ratio = uniform(rng);

      // Check ratio
      if (ratio < powerLaw(envY) / (m_linearA * envY + m_linearB)) {
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
  double powerLaw(double x) const {
    return (m_yScale * std::pow(x + m_yShift, m_yPower) + m_yPedestal) /
           m_yNorm;
  }

  Config m_cfg;

  double m_yBoundLow;
  double m_yBoundHigh;

  double m_xBoundLow;
  double m_xBoundHigh;

  double m_yScale;
  double m_yShift;
  double m_yPedestal;
  double m_yPower;
  double m_yNorm;

  std::vector<double> m_zProbs;
  std::vector<double> m_zPositions;

  std::vector<double> m_is;
  std::vector<double> m_ws;

  double m_linearA;
  double m_linearB;
};

inline E320CptBkgVertexGenerator::Config cptBkgVertexGenConfig() {
  E320Geometry::GeometryOptions gOpt;
  E320CptBkgVertexGenerator::Config cfg;

  cfg.yBoundLow = 90.3 - 29.94176 / 2;
  cfg.yBoundHigh = 331.1 + 29.94176 / 2;

  cfg.xBoundLow = gOpt.chipX - gOpt.chipSizeX / 2;
  cfg.xBoundHigh = gOpt.chipX + gOpt.chipSizeX / 2;

  cfg.yScale = 5.31456e3;
  cfg.yShift = -6.98928e1;
  cfg.yPedestal = -4.99682e3;
  cfg.yPower = -0.01;
  cfg.yNorm = 20011.3;

  cfg.zProbs = {0.29, 0.26, 0.24, 0.21};
  cfg.zPositions = {gOpt.staveZ.at(0), gOpt.staveZ.at(1), gOpt.staveZ.at(2),
                    gOpt.staveZ.at(3)};
  return cfg;
};

}  // namespace E320Sim
