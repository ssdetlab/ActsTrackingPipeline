#pragma once

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Utilities/VectorHelpers.hpp>

#include "TFile.h"
#include "TH3.h"
#include "TRandom.h"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"

/// @brief Class that samples momentum from a ROOT histogram
class HistMomentumGenerator : public IMomentumGenerator {
 public:
  struct Config {
    /// Path to file with generative
    /// histograms
    std::string pathToHist;
    /// Size histogram name
    std::string histName;
    /// Transform to global reference
    /// frame
    Acts::Transform3 transform;
  };

  HistMomentumGenerator(const Config& cfg) : m_cfg(cfg) {
    m_file = new TFile(m_cfg.pathToHist.c_str());

    m_rng = new TRandom();
    m_rng->SetSeed(std::chrono::system_clock::now().time_since_epoch().count());

    m_genPhiThetaE = (TH3D*)m_file->Get(m_cfg.histName.c_str());
  };

  ~HistMomentumGenerator() { m_file->Close(); }

  Acts::Vector3 genMomentum(RandomEngine& /*rng*/) const override {
    double phi;
    double theta;
    double E;

    m_genPhiThetaE->GetRandom3(phi, theta, E, m_rng);

    Acts::Vector3 dir{std::sin(theta) * std::cos(phi),
                      std::sin(theta) * std::sin(phi), std::cos(theta)};

    dir = m_gOpt.actsToWorld.rotation().inverse() * dir;

    return E * dir;
  }

 private:
  Config m_cfg;

  TH3D* m_genPhiThetaE = nullptr;

  TRandom* m_rng = nullptr;

  TFile* m_file = nullptr;

  E320Geometry::GeometryOptions m_gOpt;
};
