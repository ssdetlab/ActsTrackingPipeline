#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/detail/Subspace.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>

/// @brief Extended source link holding the particle hits
/// and the angles in the track coordinate system
class ExtendedSourceLink {
 public:
  /// Subspace sizes
  static constexpr const std::size_t localSubspaceSize = 4;
  static constexpr const std::size_t globalSubspaceSize = 6;

  /// @brief Surface accessor struct
  struct SurfaceAccessor {
    const Acts::Experimental::Detector* detector = nullptr;

    const Acts::Surface* operator()(const Acts::SourceLink& sourceLink) const {
      const auto& sl = sourceLink.get<ExtendedSourceLink>();
      return detector->findSurface(sl.geometryId());
    }
  };

  /// @brief Constructor
  ///
  /// @param paramsLoc particle hit and angles in the track coordinate system
  /// @param paramsGlob particle hit and momentum direction in the global frame
  /// @param cov covariance of the measurement in the track coordinate system
  /// @param gid geometry ID of the measurement surface
  /// @param eid event ID
  /// @param idx user-assigned index for fast-sim fast access
  /// @param backwards flag indicating reverse propagation of the track
  ExtendedSourceLink(const Acts::ActsVector<localSubspaceSize>& paramsLoc,
                     const Acts::ActsVector<globalSubspaceSize>& paramsGlob,
                     const Acts::ActsSquareMatrix<localSubspaceSize>& cov,
                     const Acts::GeometryIdentifier& gid, int eid, int idx,
                     bool backwards)
      : m_geometryId(gid),
        m_eventId(eid),
        m_index(idx),
        m_parametersLoc(paramsLoc),
        m_parametersGlob(paramsGlob),
        m_covariance(cov),
        m_backwards(backwards) {}

  /// Delete default-construct to satisfy SourceLinkConcept.
  ExtendedSourceLink() = delete;
  ExtendedSourceLink(const ExtendedSourceLink&) = default;
  ExtendedSourceLink(ExtendedSourceLink&&) = default;
  ExtendedSourceLink& operator=(const ExtendedSourceLink&) = default;
  ExtendedSourceLink& operator=(ExtendedSourceLink&&) = default;

  /// @brief Equality operator
  bool operator==(const ExtendedSourceLink& rhs) const {
    return (m_geometryId == rhs.geometryId()) && (m_eventId == rhs.eventId()) &&
           (m_indices == rhs.indices()) &&
           (m_parametersLoc == rhs.parametersLoc()) &&
           (m_covariance == rhs.covariance());
  }

  /// @brief Inequality operator
  bool operator!=(const ExtendedSourceLink& rhs) const {
    return !(*this == rhs);
  }

  /// @brief Get source link measurement indices
  ///
  /// @return array with indices indicating the measurement axes
  std::array<Acts::BoundIndices, localSubspaceSize> indices() const {
    return m_indices;
  }

  /// @brief Get user-assigned index
  ///
  /// @return user-assigned index
  int index() const { return m_index; }

  /// @brief Get event ID
  ///
  /// @return event ID
  int eventId() const { return m_eventId; }

  /// @brief Get geometry ID
  ///
  /// @return geometry ID
  Acts::GeometryIdentifier geometryId() const { return m_geometryId; }

  /// @brief Get source link parameters in the track frame
  ///
  /// @return source link parameters in the track frame
  Acts::ActsVector<localSubspaceSize> parametersLoc() const {
    return m_parametersLoc;
  }

  /// @brief Get source link parameters in the global frame
  ///
  /// @return source link parameters in the global frame
  Acts::ActsVector<globalSubspaceSize> parametersGlob() const {
    return m_parametersGlob;
  }

  /// @brief Get source link covariance in the track frame
  ///
  /// @return source link covariance in the track frame
  Acts::ActsSquareMatrix<localSubspaceSize> covariance() const {
    return m_covariance;
  }

  /// @brief Get source link bacwards propagation flag
  ///
  /// @return source link bacwards propagation flag
  bool isBackwards() const { return m_backwards; }

  /// @brief Set source link user-assigned index
  void setIndex(int idx) { m_index = idx; }

  /// @brief Set source link event ID
  void setEventId(int eid) { m_eventId = eid; }

 private:
  /// Geometry identifier
  Acts::GeometryIdentifier m_geometryId;

  /// Event identifier
  int m_eventId = 0u;

  /// Index for enumeration within event
  int m_index = 0u;

  /// Indices of the local coordinates
  std::array<Acts::BoundIndices, localSubspaceSize> m_indices = {
      Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundPhi, Acts::eBoundTheta};

  /// Local hit coordinates
  Acts::ActsVector<localSubspaceSize> m_parametersLoc;

  /// Global hit coordinates
  Acts::ActsVector<globalSubspaceSize> m_parametersGlob;

  /// Covariance matrix
  Acts::ActsSquareMatrix<localSubspaceSize> m_covariance;

  /// Backwards propagation flag
  bool m_backwards;
};

/// Extract the measurement from an ExtendedSourceLink.
///
/// @param gctx unused
/// @param cctx unused
/// @param sourceLink source link with the measurement
/// @param trackState TrackState to be calibrated
template <typename trajectory_t>
void extendedSourceLinkCalibratorReturn(
    const Acts::GeometryContext& /*gctx*/,
    const Acts::CalibrationContext& /*cctx*/,
    const Acts::SourceLink& sourceLink,
    typename trajectory_t::TrackStateProxy trackState) {
  const auto& sl = sourceLink.template get<ExtendedSourceLink>();

  trackState.setUncalibratedSourceLink(sourceLink);

  trackState.allocateCalibrated(ExtendedSourceLink::localSubspaceSize);
  trackState.template calibrated<ExtendedSourceLink::localSubspaceSize>() =
      sl.parametersLoc();
  trackState
      .template calibratedCovariance<ExtendedSourceLink::localSubspaceSize>() =
      sl.covariance();
  const auto& indices = sl.indices();
  trackState.setProjector(
      Acts::detail::FixedSizeSubspace<Acts::BoundIndices::eBoundSize,
                                      ExtendedSourceLink::localSubspaceSize>(
          std::array{
              indices[0],
              indices[1],
              indices[2],
              indices[3],
          })
          .projector<double>());
}

/// Extract the measurement from an ExtendedSourceLink.
///
/// @param gctx unused
/// @param cctx unused
/// @param sourceLink source link with the measurement
/// @param trackState TrackState to calibrated
template <typename trajectory_t>
void extendedSourceLinkCalibrator(
    const Acts::GeometryContext& gctx, const Acts::CalibrationContext& cctx,
    const Acts::SourceLink& sourceLink,
    typename trajectory_t::TrackStateProxy trackState) {
  extendedSourceLinkCalibratorReturn<trajectory_t>(gctx, cctx, sourceLink,
                                                   trackState);
}

/// Extract the measurement from an ExtendedSourceLink.
///
/// @param gctx unused
/// @param cctx unused
/// @param sourceLink source link with the measurement
/// @param trackState TrackState to be calibrated
template <typename trajectory_t>
void extendedSourceLinkBackwardsPhiCorrectionCalibratorReturn(
    const Acts::GeometryContext& /*gctx*/,
    const Acts::CalibrationContext& /*cctx*/,
    const Acts::SourceLink& sourceLink,
    typename trajectory_t::TrackStateProxy trackState) {
  const auto& sl = sourceLink.template get<ExtendedSourceLink>();
  Acts::ActsVector<ExtendedSourceLink::localSubspaceSize> loc =
      sl.parametersLoc();
  if (sl.isBackwards()) {
    double measPhi = loc(Acts::eBoundPhi);
    double trackPhi = trackState.predicted()(Acts::eBoundPhi);
    double phiDiff = measPhi - trackPhi;
    if (phiDiff < 0 && std::abs(phiDiff) > M_PI) {
      loc(Acts::eBoundPhi) += 2 * M_PI;
    }
    if (phiDiff > M_PI) {
      loc(Acts::eBoundPhi) -= 2 * M_PI;
    }
  }

  trackState.setUncalibratedSourceLink(sourceLink);

  trackState.allocateCalibrated(ExtendedSourceLink::localSubspaceSize);
  trackState.template calibrated<ExtendedSourceLink::localSubspaceSize>() = loc;
  trackState
      .template calibratedCovariance<ExtendedSourceLink::localSubspaceSize>() =
      sl.covariance();
  const auto& indices = sl.indices();
  trackState.setProjector(
      Acts::detail::FixedSizeSubspace<Acts::BoundIndices::eBoundSize,
                                      ExtendedSourceLink::localSubspaceSize>(
          std::array{
              indices[0],
              indices[1],
              indices[2],
              indices[3],
          })
          .projector<double>());
}

/// Extract the measurement from an ExtendedSourceLink.
///
/// @param gctx unused
/// @param cctx unused
/// @param sourceLink source link with the measurement
/// @param trackState TrackState to calibrated
template <typename trajectory_t>
void extendedSourceLinkBackwardsPhiCorrectionCalibrator(
    const Acts::GeometryContext& gctx, const Acts::CalibrationContext& cctx,
    const Acts::SourceLink& sourceLink,
    typename trajectory_t::TrackStateProxy trackState) {
  extendedSourceLinkBackwardsPhiCorrectionCalibratorReturn<trajectory_t>(
      gctx, cctx, sourceLink, trackState);
}
