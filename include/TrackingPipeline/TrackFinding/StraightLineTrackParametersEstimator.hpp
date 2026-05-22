#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/SourceLink.hpp"

#include <cstddef>

#include "TrackingPipeline/TrackFinding/ITrackParametersEstimator.hpp"
#include "TrackingPipeline/TrackFitting/StraightLineGX2Fitter.hpp"

class StraightLineTrackParametersEstimator : public ITrackParametersEstimator {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Fast GX2 fitter
    std::shared_ptr<StraightLineGX2Fitter> gx2Fitter;
    /// Number of fit iterations
    std::size_t nIterations;
    /// Maximum chi2
    double maxChi2;
    /// Reference surface
    const Acts::Surface* referenceSurface;
    /// Absolute momentum to assign
    double absMomentum;
    /// Particle hypothesis
    Acts::ParticleHypothesis particleHypothesis;
    /// Particle charge to assign
    double charge;
    /// Initial track state covariance prior
    Acts::BoundMatrix originCov;
    /// RMS of the MCS angle
    double thetaRms;
  };

  /// @brief Constructor
  explicit StraightLineTrackParametersEstimator(const Config& config);

  /// @brief Execute method
  Result estimateParameters(const Acts::GeometryContext& gctx,
                            const Acts::MagneticFieldContext& mctx,
                            const std::vector<Acts::SourceLink>& sourceLinks,
                            const std::vector<std::size_t>& sourceLinkIndices,
                            const Acts::Vector3& dir,
                            const Acts::Vector3& point) const override;

  /// @brief Readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;
};
