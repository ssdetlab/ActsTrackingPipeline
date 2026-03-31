#include "TrackingPipeline/Io/RootTrackParamsReader.hpp"

#include "Acts/EventData/ParticleHypothesis.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/PdgParticle.hpp>
#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/EventData/TrackParameters.hpp>

#include <cstddef>
#include <mutex>
#include <optional>

RootTrackParamsReader::RootTrackParamsReader(const Config& config)
    : m_cfg(config) {
  if (m_cfg.filePaths.size() == 1) {
    m_file = new TFile(m_cfg.filePaths.at(0).c_str());
    m_tree = m_file->Get<TTree>(m_cfg.treeName.c_str());
  } else {
    m_chainOwner = new TChain(m_cfg.treeName.c_str());
    // Add the files to the chain
    for (const auto& path : m_cfg.filePaths) {
      m_chainOwner->Add(path.c_str());
    }
    m_tree = dynamic_cast<TTree*>(m_chainOwner);
  }

  m_tree->SetBranchAddress("positionX", &m_params.positionX);
  m_tree->SetBranchAddress("positionY", &m_params.positionY);
  m_tree->SetBranchAddress("positionZ", &m_params.positionZ);

  m_tree->SetBranchAddress("phi", &m_params.phi);
  m_tree->SetBranchAddress("theta", &m_params.theta);

  m_tree->SetBranchAddress("qOverP", &m_params.qOverP);
  m_tree->SetBranchAddress("pdgId", &m_params.pdgId);
}

std::vector<Acts::CurvilinearTrackParameters> RootTrackParamsReader::read() {
  std::scoped_lock lock{m_readMutex};

  std::vector<Acts::CurvilinearTrackParameters> trackParams;
  trackParams.reserve(m_tree->GetEntries());
  for (std::size_t i = 0; i < m_tree->GetEntries(); i++) {
    m_tree->GetEntry(i);

    Acts::Vector4 position(m_params.positionX, m_params.positionY,
                           m_params.positionZ, 0);
    position = m_cfg.transform * position;

    Acts::Vector3 dir{std::sin(m_params.theta) * std::cos(m_params.phi),
                      std::sin(m_params.theta) * std::sin(m_params.phi),
                      std::cos(m_params.theta)};
    dir = m_cfg.transform * dir;

    trackParams.emplace_back(
        position, dir, m_params.qOverP, std::nullopt,
        Acts::ParticleHypothesis(
            static_cast<Acts::PdgParticle>(std::abs(m_params.pdgId))));
  }

  return trackParams;
}
