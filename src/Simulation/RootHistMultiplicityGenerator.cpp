#include "TrackingPipeline/Simulation/RootHistMultiplicityGenerator.hpp"

#include <chrono>
#include <cstddef>

RootHistMultiplicityGenerator::RootHistMultiplicityGenerator(const Config& cfg)
    : m_cfg(cfg) {
  m_file = new TFile(m_cfg.filePath.c_str());
  m_hist = m_file->Get<TH1I>(m_cfg.histName.c_str());

  std::size_t seed =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();
  m_rng = new TRandom(seed);
};

RootHistMultiplicityGenerator::~RootHistMultiplicityGenerator() {
  if (m_file != nullptr) {
    m_file->Close();
  }
}

std::size_t RootHistMultiplicityGenerator::genMultiplicity(
    RandomEngine& /*rng*/) const {
  return m_hist->GetRandom(m_rng);
}

double RootHistMultiplicityGenerator::getMultiplicityStdDev() const {
  return m_hist->GetStdDev();
}

double RootHistMultiplicityGenerator::getMultiplicityMean() const {
  return m_hist->GetMean();
}
