#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <random>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

/// @brief Generator reading sim particle vertex and momentum direction from a file
///
/// The generator takes a collection of files as an input and for each event
/// reads the corresponding vertex position and momentum direction from the
/// tree. The absolute momentum value is being drawn from a uniform distribution
/// on a user-defined range
class RootVertexDirectionReaderGenerator : public IVertexGenerator,
                                           public IMomentumGenerator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Input files paths
    std::vector<std::string> filePaths;
    /// Input tree name
    std::string treeName;
    /// Name of the vertex position branch
    std::string vertexBranch;
    /// Name of the momentum direction branch
    std::string directionBranch;
    /// Abolute momentum range
    double absMomentumMin;
    double absMomentumMax;
    /// Number of tree entries to skip
    std::size_t skip;
  };

  /// @brief State tracking struct
  struct State {
    /// Last generated vertex
    Acts::Vector3 vertex;
    /// Last generated momentum
    Acts::Vector3 momentum;
    /// Flag pair indicating whether the last generated
    /// vertex-momentum pair has been utilised
    std::pair<bool, bool> genState;
  };

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  explicit RootVertexDirectionReaderGenerator(const Config& cfg);

  /// @brief generate vertex
  ///
  /// @param rng random engine
  ///
  /// @return generated vertex
  Acts::Vector3 genVertex(RandomEngine& rng) const override;

  /// @brief generate momentum
  ///
  /// @param rng random engine
  ///
  /// @return generated momentum
  Acts::Vector3 genMomentum(RandomEngine& rng) const override;

  /// @brief get covariance matrix of the vertex
  Acts::SquareMatrix3 getVertexCovariance() const override;

  /// @brief get covariance matrix of the momentum
  Acts::SquareMatrix4 getMomentumCovariance() const override;

  /// @brief get vertex mean value
  Acts::Vector3 getVertexMean() const override;

  /// @brief get momentum mean value
  Acts::Vector3 getMomentumMean() const override;

 private:
  /// @brief update the internal state with a newly generated vertex-momentum pair
  ///
  /// @param rng random engine
  void internalUpdate(RandomEngine& rng) const;

  /// Configuration
  Config m_cfg;

  /// The input tree name
  TFile* m_file = nullptr;
  TTree* m_tree = nullptr;
  TChain* m_chainOwner = nullptr;

  /// Momentum
  TLorentzVector* m_originMomentum = nullptr;

  /// Vertex
  TVector3* m_vertex = nullptr;

  /// Event ID
  std::size_t m_eventId = 0;

  /// Number of tree entries
  std::size_t m_nEntries;

  /// Internal state instance
  std::unique_ptr<State> m_state;

  /// Current tree entry number
  mutable std::size_t m_idx;

  /// Uniform distribution instance
  mutable std::uniform_real_distribution<> m_uniform{0, 1};
};
