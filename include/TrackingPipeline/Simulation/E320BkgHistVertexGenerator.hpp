#pragma once

#include <string>

#include "TFile.h"
#include "TH3.h"
#include "TRandom.h"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

namespace E320Sim {

using namespace Acts::UnitLiterals;

/// @brief Class that samples vertex from a ROOT histogram
class E320BkgHistVertexGenerator : public IVertexGenerator {
 public:
  struct Config {
    /// Path to file with generative
    /// histograms
    std::string pathToHist;
    /// Size histogram name
    std::string histName;
  };

  E320BkgHistVertexGenerator(const Config& cfg) : m_cfg(cfg) {
    m_file = new TFile(m_cfg.pathToHist.c_str());

    m_rng = new TRandom();
    m_rng->SetSeed(std::chrono::system_clock::now().time_since_epoch().count());

    m_genXYZ = (TH3D*)m_file->Get(m_cfg.histName.c_str());
  };

  ~E320BkgHistVertexGenerator() { m_file->Close(); }

  Acts::Vector3 genVertex(RandomEngine& /*rng*/) const override {
    // Generate x y
    double x;
    double y;
    double z;

    m_genXYZ->GetRandom3(x, y, z, m_rng);

    // z-axis indexes layers in the histogram
    z = m_gOpt.staveZ.at(z);
    Acts::Vector3 glob(x, y, z);

    glob = m_gOpt.actsToWorld.rotation().inverse() * glob;

    return glob;
  }

 private:
  Config m_cfg;

  TH3D* m_genXYZ = nullptr;

  TRandom* m_rng = nullptr;

  TFile* m_file = nullptr;

  E320Geometry::GeometryOptions m_gOpt;
};

}  // namespace E320Sim