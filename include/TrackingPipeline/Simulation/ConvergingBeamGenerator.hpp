#pragma once

#include <Acts/Definitions/Algebra.hpp>

#include <cstddef>

#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

class ConvergingBeamGenerator : public IVertexGenerator,
                                public IMomentumGenerator {
 public:
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
    /// Angular divergence of the beam at the waist
    double waistSigmaThetaLong;
    double waistSigmaThetaShort;
    /// Beam energy
    double beamEnergy;
  };

  struct State {
    Acts::Vector3 vertex;
    Acts::Vector3 momentum;
    std::pair<bool, bool> genState;
  };

  explicit ConvergingBeamGenerator(const Config& cfg);

  Acts::Vector3 genVertex(RandomEngine& rng) const override;

  Acts::Vector3 genMomentum(RandomEngine& rng) const override;

  Acts::SquareMatrix3 getVertexCovariance() const override;

  Acts::SquareMatrix4 getMomentumCovariance() const override;

  Acts::Vector3 getVertexMean() const override;

  Acts::Vector3 getMomentumMean() const override;

 private:
  Config m_cfg;

  void internalUpdate(RandomEngine& rng) const;

  std::unique_ptr<State> m_state;

  mutable std::normal_distribution<double> m_unitNormal{0.0, 1.0};
};
