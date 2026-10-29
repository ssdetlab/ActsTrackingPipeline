#include "TrackingPipeline/Simulation/RootVertexDirectionReaderGenerator.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <cmath>
#include <stdexcept>

RootVertexDirectionReaderGenerator::RootVertexDirectionReaderGenerator(
    const Config& cfg)
    : m_cfg(cfg), m_state(std::make_unique<State>()) {
  m_state->genState = {false, false};

  if (m_cfg.filePaths.empty()) {
    throw std::invalid_argument("Missing filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

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

  //------------------------------------------------------------------
  // Set the the branches
  m_tree->SetBranchAddress(m_cfg.directionBranch.c_str(), &m_originMomentum);
  m_tree->SetBranchAddress(m_cfg.vertexBranch.c_str(), &m_vertex);
  m_nEntries = m_tree->GetEntries();
  m_idx = m_cfg.skip;
}

Acts::Vector3 RootVertexDirectionReaderGenerator::genVertex(
    RandomEngine& rng) const {
  if (!m_state->genState.first) {
    internalUpdate(rng);
  }
  m_state->genState.first = false;
  return m_state->vertex;
};

Acts::Vector3 RootVertexDirectionReaderGenerator::genMomentum(
    RandomEngine& rng) const {
  if (!m_state->genState.second) {
    internalUpdate(rng);
  }
  m_state->genState.second = false;
  return m_state->momentum;
};

Acts::SquareMatrix3 RootVertexDirectionReaderGenerator::getVertexCovariance()
    const {
  return Acts::SquareMatrix3::Identity();
}

Acts::SquareMatrix4 RootVertexDirectionReaderGenerator::getMomentumCovariance()
    const {
  return Acts::SquareMatrix4::Identity();
}

Acts::Vector3 RootVertexDirectionReaderGenerator::getVertexMean() const {
  return Acts::Vector3::Zero();
}

Acts::Vector3 RootVertexDirectionReaderGenerator::getMomentumMean() const {
  return Acts::Vector3::Zero();
}

void RootVertexDirectionReaderGenerator::internalUpdate(
    RandomEngine& rng) const {
  if (m_idx >= m_nEntries) {
    throw std::runtime_error("Reached end of the generator file");
  }
  m_tree->GetEntry(m_idx);

  m_state->vertex = Acts::Vector3(m_vertex->X(), m_vertex->Y(), m_vertex->Z());

  m_state->momentum = Acts::Vector3(
      m_originMomentum->Px(), m_originMomentum->Py(), m_originMomentum->Pz());
  m_state->momentum.normalize();

  double energy =
      m_cfg.absMomentumMin +
      (m_cfg.absMomentumMax - m_cfg.absMomentumMin) * m_uniform(rng);

  m_state->momentum *= energy;

  m_state->genState = {true, true};
  m_idx++;
}
