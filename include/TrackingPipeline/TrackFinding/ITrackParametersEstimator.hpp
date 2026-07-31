#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/Result.hpp"

/// @brief interface class for track parameters estimators
class ITrackParametersEstimator {
 public:
  /// Parameters estimation result
  using Result = Acts::Result<Acts::CurvilinearTrackParameters>;

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
  virtual Result estimateParameters(
      const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
      const std::vector<Acts::SourceLink>& sourceLinks,
      const std::vector<std::size_t>& sourceLinkIndices,
      const Acts::Vector3& dir, const Acts::Vector3& point) const = 0;
};
