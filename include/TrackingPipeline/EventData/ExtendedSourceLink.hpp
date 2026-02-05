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
#include <cstddef>

/// @brief Extended source link implementation
class ExtendedSourceLink {
 public:
  static constexpr std::size_t localSubspaceSize = 5;
  static constexpr std::size_t globalSubspaceSize = 7;

  struct SurfaceAccessor {
    const Acts::Experimental::Detector* detector = nullptr;

    const Acts::Surface* operator()(const Acts::SourceLink& sourceLink) const {
      const auto& sl = sourceLink.get<ExtendedSourceLink>();
      return detector->findSurface(sl.geometryId());
    }
  };

  /// Construct a 2d source link
  ExtendedSourceLink(const Acts::ActsVector<localSubspaceSize>& paramsLoc,
                     const Acts::ActsVector<globalSubspaceSize>& paramsGlob,
                     const Acts::ActsSquareMatrix<localSubspaceSize>& cov,
                     const Acts::GeometryIdentifier& gid, int eid, int idx)
      : m_geometryId(gid),
        m_eventId(eid),
        m_index(idx),
        m_parametersLoc(paramsLoc),
        m_parametersGlob(paramsGlob),
        m_covariance(cov) {}

  /// Default-construct an invalid source link to satisfy SourceLinkConcept.
  ExtendedSourceLink() = delete;
  ExtendedSourceLink(const ExtendedSourceLink&) = default;
  ExtendedSourceLink(ExtendedSourceLink&&) = default;
  ExtendedSourceLink& operator=(const ExtendedSourceLink&) = default;
  ExtendedSourceLink& operator=(ExtendedSourceLink&&) = default;

  bool operator==(const ExtendedSourceLink& rhs) const {
    return (m_geometryId == rhs.geometryId()) && (m_eventId == rhs.eventId()) &&
           (m_indices == rhs.indices()) &&
           (m_parametersLoc == rhs.parametersLoc()) &&
           (m_covariance == rhs.covariance());
  }

  bool operator!=(const ExtendedSourceLink& rhs) const {
    return !(*this == rhs);
  }

  std::array<Acts::BoundIndices, localSubspaceSize> indices() const {
    return m_indices;
  }

  int index() const { return m_index; }

  int eventId() const { return m_eventId; }

  Acts::GeometryIdentifier geometryId() const { return m_geometryId; }

  Acts::ActsVector<localSubspaceSize> parametersLoc() const {
    return m_parametersLoc;
  }

  Acts::ActsVector<globalSubspaceSize> parametersGlob() const {
    return m_parametersGlob;
  }

  Acts::ActsSquareMatrix<localSubspaceSize> covariance() const {
    return m_covariance;
  }

  void setIndex(int idx) { m_index = idx; }

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
      Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundPhi, Acts::eBoundTheta,
      Acts::eBoundQOverP};

  /// Local hit coordinates
  Acts::ActsVector<localSubspaceSize> m_parametersLoc;

  /// Global hit coordinates
  Acts::ActsVector<globalSubspaceSize> m_parametersGlob;

  /// Covariance matrix
  Acts::ActsSquareMatrix<localSubspaceSize> m_covariance;
};

/// Extract the measurement from a ExtendedSourceLink.
///
/// @param gctx Unused
/// @param trackState TrackState to calibrated
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
          std::array{indices[0], indices[1], indices[2], indices[3],
                     indices[4]})
          .projector<double>());
}

/// Extract the measurement from a ExtendedSourceLink.
///
/// @param gctx Unused
/// @param trackState TrackState to calibrated
template <typename trajectory_t>
void extendedSourceLinkCalibrator(
    const Acts::GeometryContext& gctx, const Acts::CalibrationContext& cctx,
    const Acts::SourceLink& sourceLink,
    typename trajectory_t::TrackStateProxy trackState) {
  extendedSourceLinkCalibratorReturn<trajectory_t>(gctx, cctx, sourceLink,
                                                   trackState);
}
