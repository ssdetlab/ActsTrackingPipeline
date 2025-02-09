#pragma once

#include "Acts/Definitions/Units.hpp"

#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TRandom.h"
#include "TrackingPipeline/Simulation/IDigitizer.hpp"

namespace E320Sim {
using namespace Acts::UnitLiterals;

/// @brief Class that digitizes hits based on the
/// analysis of the E320 experiment data
class E320HistDigitizer : public IDigitizer {
 public:
  struct Config {
    /// Path to file with generative
    /// histograms
    std::string pathToHist;
    /// Size histogram name
    std::string histName;
  };

  E320HistDigitizer(const Config& cfg) : m_cfg(cfg) {
    m_file = new TFile(m_cfg.pathToHist.c_str());

    m_rng = new TRandom();
    m_rng->SetSeed(std::chrono::system_clock::now().time_since_epoch().count());

    m_genSize = (TH1D*)m_file->Get(m_cfg.histName.c_str());
  };

  ~E320HistDigitizer() { m_file->Close(); }

  std::pair<Acts::SquareMatrix2, Acts::Vector2> genCluster(
      RandomEngine& /*rng*/, Acts::GeometryIdentifier /*geoId*/,
      Acts::Vector2 pos) const override {
    int size = m_genSize->GetRandom(m_rng);

    Acts::ActsScalar errX = m_pixSizeX / std::sqrt(12 * size);
    Acts::ActsScalar errY = m_pixSizeY / std::sqrt(12 * size);
    Acts::Vector2 stdDev(errX, errY);
    Acts::SquareMatrix2 cov = stdDev.cwiseProduct(stdDev).asDiagonal();

    return {cov, pos};
  }

 private:
  Config m_cfg;

  Acts::ActsScalar m_pixSizeX = 27_um;
  Acts::ActsScalar m_pixSizeY = 29_um;

  TH1D* m_genSize = nullptr;

  TRandom* m_rng = nullptr;

  TFile* m_file = nullptr;
};

}  // namespace E320Sim
