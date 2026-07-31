#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <random>

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

/// @brief Vertex and momentum generator simulating a converging beam
///
/// The generator simulates a beam originating at a specified reference surface
/// and passing through a specified waist plane with a given shot-to-shot
/// spatial and angular jitter and a unifomly distributed in a specified range
/// energy
class ConvergingBeamGenerator : public IVertexGenerator,
                                public IMomentumGenerator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Indices
    std::size_t primaryIdx;
    std::size_t longIdx;
    std::size_t shortIdx;
    /// Reference plane
    double referencePositionPrimary;
    /// Beam waist position
    Acts::Vector3 waistPosition;
    /// Beam waist size
    double waistSigmaLong;
    double waistSigmaShort;
    /// Angle of the beam at the waist
    double waistMeanThetaLong;
    double waistMeanThetaShort;
    /// Angular divergence of the beam at the waist
    double waistSigmaThetaLong;
    double waistSigmaThetaShort;
    /// Beam energy range
    double beamEnergyMin;
    double beamEnergyMax;
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
  explicit ConvergingBeamGenerator(const Config& cfg);

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

  /// Internal state instance
  std::unique_ptr<State> m_state;

  /// Normal distribution instance
  mutable std::normal_distribution<double> m_normal{0.0, 1.0};

  /// Uniform distribution instance
  mutable std::uniform_real_distribution<double> m_uniform{0.0, 1.0};
};
