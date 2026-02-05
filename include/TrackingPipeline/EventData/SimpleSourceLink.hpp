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

/// @brief A minimal source link implementation
/// that is easy to convert to the Measurement
///
/// @note Stores the geometry identifier,
/// local hit coordinates, and the covariance
class SimpleSourceLink {
 public:
  static constexpr std::size_t localSubspaceSize = 2;
  static constexpr std::size_t globalSubspaceSize = 3;

  struct SurfaceAccessor {
    const Acts::Experimental::Detector* detector = nullptr;

    const Acts::Surface* operator()(const Acts::SourceLink& sourceLink) const {
      const auto& sl = sourceLink.get<SimpleSourceLink>();
      return detector->findSurface(sl.geometryId());
    }
  };

  /// Construct a 2d source link
  SimpleSourceLink(const Acts::ActsVector<localSubspaceSize>& paramsLoc,
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
  SimpleSourceLink() = delete;
  SimpleSourceLink(const SimpleSourceLink&) = default;
  SimpleSourceLink(SimpleSourceLink&&) = default;
  SimpleSourceLink& operator=(const SimpleSourceLink&) = default;
  SimpleSourceLink& operator=(SimpleSourceLink&&) = default;

  bool operator==(const SimpleSourceLink& rhs) const {
    return (m_geometryId == rhs.geometryId()) && (m_eventId == rhs.eventId()) &&
           (m_indices == rhs.indices()) &&
           (m_parametersLoc == rhs.parametersLoc()) &&
           (m_covariance == rhs.covariance());
  }

  bool operator!=(const SimpleSourceLink& rhs) const { return !(*this == rhs); }

  std::array<Acts::BoundIndices, localSubspaceSize> indices() const {
    return m_indices;
  }

  int index() const { return m_index; }

  int eventId() const { return m_eventId; }

  Acts::GeometryIdentifier geometryId() const { return m_geometryId; }

  Acts::Vector2 parametersLoc() const { return m_parametersLoc; }

  Acts::Vector3 parametersGlob() const { return m_parametersGlob; }

  Acts::SquareMatrix2 covariance() const { return m_covariance; }

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
      Acts::eBoundLoc0, Acts::eBoundLoc1};

  /// Local hit coordinates
  Acts::ActsVector<localSubspaceSize> m_parametersLoc;

  /// Global hit coordinates
  Acts::ActsVector<globalSubspaceSize> m_parametersGlob;

  /// Covariance matrix
  Acts::ActsSquareMatrix<localSubspaceSize> m_covariance;
};

/// Extract the measurement from a SimpleSourceLink.
///
/// @param gctx Unused
/// @param trackState TrackState to calibrated
template <typename trajectory_t>
void simpleSourceLinkCalibratorReturn(
    const Acts::GeometryContext& /*gctx*/,
    const Acts::CalibrationContext& /*cctx*/,
    const Acts::SourceLink& sourceLink,
    typename trajectory_t::TrackStateProxy trackState) {
  const auto& sl = sourceLink.template get<SimpleSourceLink>();

  trackState.setUncalibratedSourceLink(sourceLink);

  trackState.allocateCalibrated(SimpleSourceLink::localSubspaceSize);
  trackState.template calibrated<SimpleSourceLink::localSubspaceSize>() =
      sl.parametersLoc();
  trackState
      .template calibratedCovariance<SimpleSourceLink::localSubspaceSize>() =
      sl.covariance();
  const auto& indices = sl.indices();
  trackState.setProjector(
      Acts::detail::FixedSizeSubspace<Acts::BoundIndices::eBoundSize,
                                      SimpleSourceLink::localSubspaceSize>(
          std::array{indices[0], indices[1]})
          .projector<double>());
}

/// Extract the measurement from a SimpleSourceLink.
///
/// @param gctx Unused
/// @param trackState TrackState to calibrated
template <typename trajectory_t>
void simpleSourceLinkCalibrator(
    const Acts::GeometryContext& gctx, const Acts::CalibrationContext& cctx,
    const Acts::SourceLink& sourceLink,
    typename trajectory_t::TrackStateProxy trackState) {
  simpleSourceLinkCalibratorReturn<trajectory_t>(gctx, cctx, sourceLink,
                                                 trackState);
}
