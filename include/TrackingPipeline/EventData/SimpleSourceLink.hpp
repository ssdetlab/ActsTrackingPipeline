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

/// @brief Simple source link holding the particle hits
/// in the track coordinate system
class SimpleSourceLink {
 public:
  /// Subspace sizes
  static constexpr const std::size_t localSubspaceSize = 2;
  static constexpr const std::size_t globalSubspaceSize = 3;

  /// @brief Surface accessor struct
  struct SurfaceAccessor {
    const Acts::Experimental::Detector* detector = nullptr;

    const Acts::Surface* operator()(const Acts::SourceLink& sourceLink) const {
      const auto& sl = sourceLink.get<SimpleSourceLink>();
      return detector->findSurface(sl.geometryId());
    }
  };

  /// @brief Constructor
  ///
  /// @param paramsLoc particle hit in the track coordinate system
  /// @param paramsGlob particle hit in the global frame
  /// @param cov covariance of the measurement in the track coordinate system
  /// @param gid geometry ID of the measurement surface
  /// @param eid event ID
  /// @param idx user-assigned index for fast-sim fast access
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

  /// Delete default-construct to satisfy SourceLinkConcept.
  SimpleSourceLink() = delete;
  SimpleSourceLink(const SimpleSourceLink&) = default;
  SimpleSourceLink(SimpleSourceLink&&) = default;
  SimpleSourceLink& operator=(const SimpleSourceLink&) = default;
  SimpleSourceLink& operator=(SimpleSourceLink&&) = default;

  /// @brief Equality operator
  bool operator==(const SimpleSourceLink& rhs) const {
    return (m_geometryId == rhs.geometryId()) && (m_eventId == rhs.eventId()) &&
           (m_indices == rhs.indices()) &&
           (m_parametersLoc == rhs.parametersLoc()) &&
           (m_covariance == rhs.covariance());
  }

  /// @brief Inequality operator
  bool operator!=(const SimpleSourceLink& rhs) const { return !(*this == rhs); }

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
  Acts::Vector2 parametersLoc() const { return m_parametersLoc; }

  /// @brief Get source link parameters in the global frame
  ///
  /// @return source link parameters in the global frame
  Acts::Vector3 parametersGlob() const { return m_parametersGlob; }

  /// @brief Get source link covariance in the track frame
  ///
  /// @return source link covariance in the track frame
  Acts::SquareMatrix2 covariance() const {
    return m_covarianceAnnealingFactor * m_covariance;
  }

  /// @brief Set source link user-assigned index
  void setIndex(int idx) { m_index = idx; }

  /// @brief Set source link event ID
  void setEventId(int eid) { m_eventId = eid; }

  /// @brief Set annealing factor for the covariance
  void setCovarianceAnnealingFactor(double alpha) {
    m_covarianceAnnealingFactor = alpha;
  }

  /// @brief Set source link covariance
  void setCovariance(const Acts::SquareMatrix2& cov) { m_covariance = cov; }

 private:
  /// Geometry identifier
  Acts::GeometryIdentifier m_geometryId;

  /// Event identifier
  int m_eventId = 0u;

  /// Index for enumeration within event
  int m_index = 0u;

  /// Anneling factor of the covariance matrix
  double m_covarianceAnnealingFactor = 1;

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
/// @param gctx unused
/// @param cctx unused
/// @param sourceLink source link with the measurement
/// @param trackState TrackState to be calibrated
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
/// @param gctx unused
/// @param cctx unused
/// @param sourceLink source link with the measurement
/// @param trackState TrackState to be calibrated
template <typename trajectory_t>
void simpleSourceLinkCalibrator(
    const Acts::GeometryContext& gctx, const Acts::CalibrationContext& cctx,
    const Acts::SourceLink& sourceLink,
    typename trajectory_t::TrackStateProxy trackState) {
  simpleSourceLinkCalibratorReturn<trajectory_t>(gctx, cctx, sourceLink,
                                                 trackState);
}
