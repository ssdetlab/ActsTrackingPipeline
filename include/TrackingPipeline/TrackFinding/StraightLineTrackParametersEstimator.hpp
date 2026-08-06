#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"

#include <cstddef>

#include "TrackingPipeline/TrackFinding/ITrackParametersEstimator.hpp"
#include "TrackingPipeline/TrackFitting/StraightLineGX2Fitter.hpp"

/// @brief Track parameters estimator performing a straing line track estimation
///
/// The estimator performs a global chi2 fit of the measurements to infer the
/// position and the direction of the track at the reference surface. The
/// momentum assigned to the track is specified by the user
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
  ///
  /// @param config configuration struct
  explicit StraightLineTrackParametersEstimator(const Config& config);

  /// @brief perform track parameters estimation
  ///
  /// @param gctx current geometry context
  /// @param mctx current magnetic field context
  /// @param sourceLinksContainer event source links collection
  /// @param dir initial guess for the track direction
  /// @param point initial guess for the track passing point
  ///
  /// @return esimtated track parameters in a Result wrapper
  Result estimateParameters(const Acts::GeometryContext& gctx,
                            const Acts::MagneticFieldContext& mctx,
                            const SourceLinkContainer& sourceLinkContainer,
                            const Acts::Vector3& dir,
                            const Acts::Vector3& point) const override;

  /// @brief Readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;
};
