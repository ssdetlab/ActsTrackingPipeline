#pragma once

#include "Acts/Utilities/detail/Subspace.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <stdexcept>

/// @brief A minimal source link implementation
/// that is easy to convert to the Measurement
/// 
/// @note Stores the geometry identifier,
/// local hit coordinates, and the covariance
class SimpleSourceLink {
    public: 
        struct SurfaceAccessor {
            const Acts::Experimental::Detector* detector = nullptr;
        
            const Acts::Surface* operator()(
                const Acts::SourceLink& sourceLink) const {
                    const auto& sl = sourceLink.get<SimpleSourceLink>();
                    return detector->findSurface(sl.geometryId());
            }
        };

        /// Construct a 2d source link
        SimpleSourceLink(
            const Acts::ActsVector<2>& params,
            const Acts::ActsSquareMatrix<2>& cov,
            Acts::GeometryIdentifier gid,
            std::int32_t eid,
            std::int32_t idx)
            : m_geometryId(gid),
            m_eventId(eid),
            m_index(idx),
            m_parameters(params),
            m_covariance(cov) {}
    
        /// Default-construct an invalid source link to satisfy SourceLinkConcept.
        SimpleSourceLink() = delete;
        SimpleSourceLink(const SimpleSourceLink&) = default;
        SimpleSourceLink(SimpleSourceLink&&) = default;
        SimpleSourceLink& operator=(const SimpleSourceLink&) = default;
        SimpleSourceLink& operator=(SimpleSourceLink&&) = default;

        bool operator==(const SimpleSourceLink& rhs) const {
            return (m_geometryId == rhs.geometryId()) && 
                (m_eventId == rhs.eventId()) &&
                (m_indices == rhs.indices()) && 
                (m_parameters == rhs.parameters()) &&
                (m_covariance == rhs.covariance());
        }

        bool operator!=(const SimpleSourceLink& rhs) const { return !(*this == rhs); }

        std::array<Acts::BoundIndices, 2> indices() const { return m_indices; }

        std::int32_t index() const { return m_index; }

        std::int32_t eventId() const { return m_eventId; }

        Acts::GeometryIdentifier geometryId() const { return m_geometryId; }

        Acts::Vector2 parameters() const { return m_parameters; }

        Acts::SquareMatrix2 covariance() const { return m_covariance; }

        void setIndex(std::int32_t idx) { m_index = idx; }

        void setEventId(std::int32_t eid) { m_eventId = eid; }

    private:
        /// Geometry identifier
        Acts::GeometryIdentifier m_geometryId;

        /// Event identifier
        std::int32_t m_eventId = 0u;
    
        /// Index for global matching   
        std::int32_t m_index = 0u;
    
        /// Indices of the local coordinates
        std::array<Acts::BoundIndices, 2> m_indices = 
            {Acts::eBoundLoc0, Acts::eBoundLoc1};
    
        /// Local hit coordinates
        Acts::ActsVector<2> m_parameters;
    
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
        trackState.template calibrated<2>() = sl.parameters();
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
    const Acts::GeometryContext& gctx, 
    const Acts::CalibrationContext& cctx,
    const Acts::SourceLink& sourceLink,
    typename trajectory_t::TrackStateProxy trackState) {
        simpleSourceLinkCalibratorReturn<trajectory_t>(
            gctx, cctx, sourceLink, trackState);
}

// Calibrator to transform the source links
// to global coordinates
class SimpleSourceLinkCoordinateCalibrator {
    public:
        Acts::SourceLinkSurfaceAccessor m_surfaceAccessor;

        Acts::Vector3 operator()(
            const Acts::GeometryContext& geoCtx,
            const Acts::SourceLink& sourceLink) const {
                auto ssl = sourceLink.get<SimpleSourceLink>();
                auto res = m_surfaceAccessor(sourceLink)->localToGlobal(
                    geoCtx, ssl.parameters(), Acts::Vector3{0, 1, 0});
                return res;
        }
};
