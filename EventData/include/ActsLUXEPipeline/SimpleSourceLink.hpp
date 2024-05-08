#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/EventData/Measurement.hpp"
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
struct SimpleSourceLink final {
    Acts::GeometryIdentifier m_geometryId{};
    std::int32_t eventId = 0u;
    std::array<Acts::BoundIndices, 2> indices = 
        {Acts::eBoundLoc0, Acts::eBoundLoc1};
    Acts::ActsVector<2> parameters;
    Acts::ActsSquareMatrix<2> covariance;

    /// Construct a 2d source link
    SimpleSourceLink(const Acts::ActsVector<2>& params,
        const Acts::ActsSquareMatrix<2>& cov,
        Acts::GeometryIdentifier gid,
        std::int32_t eid)
        : m_geometryId(gid),
        eventId(eid),
        parameters(params),
        covariance(cov) {}

    /// Default-construct an invalid source link to satisfy SourceLinkConcept.
    SimpleSourceLink() = default;
    SimpleSourceLink(const SimpleSourceLink&) = default;
    SimpleSourceLink(SimpleSourceLink&&) = default;
    SimpleSourceLink& operator=(const SimpleSourceLink&) = default;
    SimpleSourceLink& operator=(SimpleSourceLink&&) = default;
    
    bool operator==(const SimpleSourceLink& rhs) const {
        return (m_geometryId == rhs.m_geometryId) && (eventId == rhs.eventId) &&
            (indices == rhs.indices) && (parameters == rhs.parameters) &&
            (covariance == rhs.covariance);
    }
    bool operator!=(const SimpleSourceLink& rhs) const { return !(*this == rhs); }
    std::ostream& print(std::ostream& os) const {
        os << "SimpleSourceLink(geometryId=" << m_geometryId
            << ",eventId=" << eventId;
            os << ")";
        return os;
    }
    constexpr std::int32_t index() const { return eventId; }

    Acts::GeometryIdentifier geometryId() const { return m_geometryId; }

    struct SurfaceAccessor {
        const Acts::Experimental::Detector& detector;
    
        const Acts::Surface* operator()(const Acts::SourceLink& sourceLink) const {
            const auto& sl = sourceLink.get<SimpleSourceLink>();
            return *detector.sensitiveHierarchyMap().find(
                sl.geometryId());
        }
    };

};

inline std::ostream& operator<<(std::ostream& os,
    const SimpleSourceLink& sourceLink) {
        return sourceLink.print(os);
}

/// Extract the measurement from a SimpleSourceLink.
///
/// @param gctx Unused
/// @param trackState TrackState to calibrated
/// @return The measurement used
template <typename trajectory_t>
Acts::BoundVariantMeasurement simpleSourceLinkCalibratorReturn(
    const Acts::GeometryContext& /*gctx*/, 
    const Acts::CalibrationContext& /*cctx*/,
    const Acts::SourceLink& sourceLink,
    typename trajectory_t::TrackStateProxy trackState) {
        SimpleSourceLink sl = sourceLink.template get<SimpleSourceLink>();

        trackState.setUncalibratedSourceLink(sourceLink);

        auto meas =
            makeMeasurement(trackState.getUncalibratedSourceLink(), 
                sl.parameters, sl.covariance, 
                sl.indices[0], sl.indices[1]);
        trackState.allocateCalibrated(2);
        trackState.setCalibrated(meas);
        return meas;
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
