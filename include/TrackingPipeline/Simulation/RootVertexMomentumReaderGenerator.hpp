#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <random>
#include <string>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>

#include "TTree.h"
#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

class RootVertexMomentumReaderGenerator : public IVertexGenerator,
                                           public IMomentumGenerator {
 public:
  struct Config {
    std::vector<std::string> filePaths;
    std::string treeName;
    std::string vertexBranch;
    std::string directionBranch;
    double absMomentumMin;
    double absMomentumMax;
    std::size_t startIdx;
  };

  struct State {
    Acts::Vector3 vertex;
    Acts::Vector3 momentum;
    std::pair<bool, bool> genState;
  };

  explicit RootVertexMomentumReaderGenerator(const Config& cfg);

  Acts::Vector3 genVertex(RandomEngine& rng) const override;

  Acts::Vector3 genMomentum(RandomEngine& rng) const override;

  Acts::SquareMatrix3 getVertexCovariance() const override;

  Acts::SquareMatrix4 getMomentumCovariance() const override;

  Acts::Vector3 getVertexMean() const override;

  Acts::Vector3 getMomentumMean() const override;

 private:
  Config m_cfg;

  void internalUpdate(RandomEngine& rng) const;

  /// The input tree name
  TFile* m_file = nullptr;
  TTree* m_tree = nullptr;
  TChain* m_chainOwner = nullptr;

  /// Momentum
  TLorentzVector* m_originMomentum = nullptr;

  /// Vertex
  TVector3* m_vertex = nullptr;

  std::size_t m_eventId = 0;

  std::unique_ptr<State> m_state;

  std::size_t m_nEntries;
  mutable std::size_t m_idx;
  mutable std::uniform_real_distribution<> m_uniform{0, 1};
};
