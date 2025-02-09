#include "TrackingPipeline/Io/RootTrackParamsReader.hpp"

#include "Acts/EventData/ParticleHypothesis.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/PdgParticle.hpp>
#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/EventData/TrackParameters.hpp>

#include <cstddef>
#include <optional>

RootTrackParamsReader::RootTrackParamsReader(const Config& config)
    : m_cfg(config) {
  m_chain = new TChain(m_cfg.treeName.c_str());
  for (const auto& path : m_cfg.filePaths) {
    m_chain->Add(path.c_str());
  }

  m_chain->SetBranchAddress("positionX", &m_params.positionX);
  m_chain->SetBranchAddress("positionY", &m_params.positionY);
  m_chain->SetBranchAddress("positionZ", &m_params.positionZ);

  m_chain->SetBranchAddress("phi", &m_params.phi);
  m_chain->SetBranchAddress("theta", &m_params.theta);

  m_chain->SetBranchAddress("qOverP", &m_params.qOverP);
  m_chain->SetBranchAddress("pdgId", &m_params.pdgId);
}

std::vector<Acts::CurvilinearTrackParameters> RootTrackParamsReader::read() {
  std::lock_guard<std::mutex> lock(m_readMutex);

  std::vector<Acts::CurvilinearTrackParameters> trackParams;
  trackParams.reserve(m_chain->GetEntries());
  for (std::size_t i = 0; i < m_chain->GetEntries(); i++) {
    m_chain->GetEntry(i);

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
