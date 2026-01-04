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
  struct SurfaceAccessor {
    const Acts::Experimental::Detector* detector = nullptr;

    const Acts::Surface* operator()(const Acts::SourceLink& sourceLink) const {
      const auto& sl = sourceLink.get<SimpleSourceLink>();
      return detector->findSurface(sl.geometryId());
    }
  };

  /// Construct a 2d source link
  SimpleSourceLink(const Acts::ActsVector<2>& paramsLoc,
                   const Acts::ActsVector<3>& paramsGlob,
                   const Acts::ActsSquareMatrix<2>& cov,
                   Acts::GeometryIdentifier gid, int eid, int idx)
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

  std::array<Acts::BoundIndices, 2> indices() const { return m_indices; }

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
  std::array<Acts::BoundIndices, 2> m_indices = {Acts::eBoundLoc0,
                                                 Acts::eBoundLoc1};

  /// Local hit coordinates
  Acts::ActsVector<2> m_parametersLoc;

  /// Global hit coordinates
  Acts::ActsVector<3> m_parametersGlob;

  /// Covariance matrix
  Acts::ActsSquareMatrix<2> m_covariance;
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
  SimpleSourceLink sl = sourceLink.template get<SimpleSourceLink>();

  trackState.setUncalibratedSourceLink(sourceLink);

  trackState.allocateCalibrated(2);
  trackState.template calibrated<2>() = sl.parametersLoc();
  trackState.template calibratedCovariance<2>() = sl.covariance();
  trackState.setProjector(
      Acts::detail::FixedSizeSubspace<Acts::BoundIndices::eBoundSize, 2>(
          std::array{sl.indices()[0], sl.indices()[1]})
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

// Calibrator to transform the source links
// to global coordinates
class SimpleSourceLinkCoordinateCalibrator {
 public:
  Acts::SourceLinkSurfaceAccessor m_surfaceAccessor;

  Acts::Vector3 operator()(const Acts::GeometryContext& geoCtx,
                           const Acts::SourceLink& sourceLink) const {
    auto ssl = sourceLink.get<SimpleSourceLink>();
    auto res = m_surfaceAccessor(sourceLink)
                   ->localToGlobal(geoCtx, ssl.parametersLoc(),
                                   Acts::Vector3{0, 1, 0});
    return res;
  }
};
