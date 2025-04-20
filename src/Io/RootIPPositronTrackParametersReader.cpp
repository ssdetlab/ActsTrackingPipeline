#include "TrackingPipeline/Io/RootIPPositronTrackParametersReader.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"

#include <cstddef>
#include <mutex>
#include <optional>

using namespace Acts::UnitLiterals;

RootIPPositronTrackParametersReader::RootIPPositronTrackParametersReader(
    const Config& config)
    : m_cfg(config) {
  m_chain = new TChain(m_cfg.treeName.c_str());
  for (const auto& path : m_cfg.filePaths) {
    m_chain->Add(path.c_str());
  }

  m_chain->SetBranchAddress("vx", &m_positionX);
  m_chain->SetBranchAddress("vy", &m_positionY);
  m_chain->SetBranchAddress("vz", &m_positionZ);

  m_chain->SetBranchAddress("px", &m_phi);
  m_chain->SetBranchAddress("py", &m_theta);
  m_chain->SetBranchAddress("pz", &m_E);

  m_chain->SetBranchAddress("pdgId", &m_pdgId);
}

std::vector<Acts::CurvilinearTrackParameters>
RootIPPositronTrackParametersReader::read() {
  std::scoped_lock lock{m_readMutex};

  std::vector<Acts::CurvilinearTrackParameters> trackParams;
  trackParams.reserve(m_chain->GetEntries());
  for (std::size_t i = 0; i < m_chain->GetEntries(); i++) {
    m_chain->GetEntry(i);

    for (std::size_t j = 0; j < m_E->size(); j++) {
       Acts::Vector4 position(m_positionY->at(j), m_positionZ->at(j),
                              m_positionX->at(j), 0);
      position = Acts::Vector4(position.x(), position.z(), -position.y(), 0);
      Acts::Vector3 mom(m_phi->at(j), m_E->at(j), -m_theta->at(j));

      double qOverP = 1_e / mom.norm();
      Acts::Vector3 dir = mom / mom.norm();
      trackParams.emplace_back(
          position, dir, qOverP, std::nullopt,
          Acts::ParticleHypothesis(
              static_cast<Acts::PdgParticle>(std::abs(m_pdgId->at(j)))));
    }
  }

  return trackParams;
}
