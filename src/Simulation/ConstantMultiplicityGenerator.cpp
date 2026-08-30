#include "TrackingPipeline/Simulation/ConstantMultiplicityGenerator.hpp"

ConstantMultiplicityGenerator::ConstantMultiplicityGenerator(const Config& cfg)
    : m_cfg(cfg) {};

std::size_t ConstantMultiplicityGenerator::genMultiplicity(
    RandomEngine& /*rng*/) const {
  return m_cfg.eventMultiplicity;
}

double ConstantMultiplicityGenerator::getMultiplicityStdDev() const {
  return 0;
}

double ConstantMultiplicityGenerator::getMultiplicityMean() const {
  return m_cfg.eventMultiplicity;
}
