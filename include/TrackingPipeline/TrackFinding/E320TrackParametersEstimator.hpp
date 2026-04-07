#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Utilities/Result.hpp>

#include <cstddef>

#include "TrackingPipeline/TrackFinding/ITrackParametersEstimator.hpp"
#include "TrackingPipeline/TrackFitting/FastGX2Fitter.hpp"

namespace E320 {

class E320TrackParametersEstimator : public ITrackParametersEstimator {
 public:
  enum struct PropagationDirection : int { forward = 0, backward = 1 };

  /// @brief The nested configuration struct
  struct Config {
    /// Fast GX2 fitter
    std::shared_ptr<FastGX2Fitter> gx2Fitter;
    /// Number of fit iterations
    std::size_t nIterations;
    /// Maximum chi2
    double maxChi2;
    /// Reference surface
    const Acts::Surface* referenceSurface;
    /// Initial track state covariance prior
    Acts::BoundMatrix originCov;
    /// Forward or backward propagation parameters
    PropagationDirection propDirection;
  };

  /// @brief Constructor
  explicit E320TrackParametersEstimator(const Config& config);

  /// @brief Execute method
  Result estimateParameters(const Acts::GeometryContext& gctx,
                            const std::vector<Acts::SourceLink>& sourceLinks,
                            const std::vector<std::size_t>& sourceLinkIndices,
                            const Acts::Vector3& dir,
                            const Acts::Vector3& point) const override;

  /// @brief Readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;

  Acts::BoundMatrix transportCovToReference(
      const Acts::GeometryContext& gctx, const Acts::Vector3& refSurfacePoint,
      const Acts::Vector3& point, const Acts::Vector3& dir,
      const Acts::ActsSquareMatrix<6>& cov) const;

  double m_dipoleLength;
  double m_dipoleFieldStrength;
};

}  // namespace E320
