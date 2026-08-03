#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"

#include <cstddef>

#include "TrackingPipeline/TrackFinding/ITrackParametersEstimator.hpp"
#include "TrackingPipeline/TrackFitting/StraightLineGX2Fitter.hpp"

namespace E320 {

/// @brief e320-specific track parameters esimator
///
/// Track parameters estimator accounting for e320 specifics. The position and
/// direction estimation is done by performing a global chi2 fit of the tracker
/// hits. The fitted parameters are propagated to a specified reference surface.
/// The q/P estimation is done through extracting the dipole strength
/// information from the E320GeometryOptions. The MCS is accounted for by
/// extracting the pixel thickness from the E320GeometryOptions.
class E320TrackParametersEstimator : public ITrackParametersEstimator {
 public:
  /// Flag for foward/backward propagation parameters estimation
  enum struct PropagationDirection : int { forward = 0, backward = 1 };

  /// @brief The nested configuration struct
  struct Config {
    /// GX2 fitter
    std::shared_ptr<StraightLineGX2Fitter> gx2Fitter;
    /// Number of GX2 fit iterations
    std::size_t nIterations;
    /// Maximum GX2 chi2 cut
    double maxChi2;
    /// Reference surface
    const Acts::Surface* referenceSurface;
    /// Initial track state covariance prior
    Acts::BoundMatrix originCov;
    /// Forward or backward propagation parameters
    PropagationDirection propDirection;
  };

  /// @brief Constructor
  ///
  /// @param config configuration struct
  explicit E320TrackParametersEstimator(const Config& config);

  /// @brief perform track parameters estimation
  ///
  /// @param gctx current geometry context
  /// @param mctx current magnetic field context
  /// @param sourceLinks event source links collection
  /// @param sourceLinkIndices event source link indices collection
  /// @param dir initial guess for the track direction
  /// @param point initial guess for the track passing point
  ///
  /// @return esimtated track parameters in a Result wrapper
  Result estimateParameters(const Acts::GeometryContext& gctx,
                            const Acts::MagneticFieldContext& mctx,
                            const std::vector<Acts::SourceLink>& sourceLinks,
                            const std::vector<std::size_t>& sourceLinkIndices,
                            const Acts::Vector3& dir,
                            const Acts::Vector3& point) const override;

  /// @brief Readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// @brief transport the GX2 estimated covariance to the reference surface
  ///
  /// @param gctx current geometry context
  /// @param refSurfacePoint point on the reference surface
  /// @param point GX2 track passing point estimate
  /// @param dir GX2 track direction estimate
  /// @param cov GX2 track covariance esimtate
  /// @param dipoleFieldStrength dipole field in the current event
  ///
  /// @return track covariance estimate at the reference surface
  Acts::BoundMatrix transportCovToReference(
      const Acts::GeometryContext& gctx, const Acts::Vector3& refSurfacePoint,
      const Acts::Vector3& point, const Acts::Vector3& dir,
      const StraightLineGX2Fitter::Covariance& cov,
      double dipoleFieldStrength) const;

  /// Configuration
  Config m_cfg;

  /// Dipole field direction index
  std::size_t m_dipoleDirIdx;

  /// Dipole field ID
  std::size_t m_dipoleId;

  /// Length of the dipole
  double m_dipoleLength;

  /// Dipole field strength
  double m_dipoleFieldStrength;

  /// Pixel sensor thickness
  double m_sensorThickness;

  /// Indident particle charge
  double m_particleCharge;
};

}  // namespace E320
