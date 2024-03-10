#pragma once

#include "Acts/Detector/Detector.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <cassert>
#include <boost/container/flat_set.hpp>

struct GeometryIdGetter {
    // explicit geometry identifier are just forwarded
    constexpr Acts::GeometryIdentifier operator()(
        Acts::GeometryIdentifier geometryId) const {
        return geometryId;
    }
    // encoded geometry ids are converted back to geometry identifiers.
    constexpr Acts::GeometryIdentifier operator()(
        Acts::GeometryIdentifier::Value encoded) const {
        return Acts::GeometryIdentifier(encoded);
    }
    // support elements in map-like structures.
    template <typename T>
    constexpr Acts::GeometryIdentifier operator()(
        const std::pair<Acts::GeometryIdentifier, T>& mapItem) const {
        return mapItem.first;
    }
    // support elements that implement `.geometryId()`.
    template <typename T>
    inline auto operator()(const T& thing) const
        -> decltype(thing.geometryId(), Acts::GeometryIdentifier()) {
        return thing.geometryId();
    }
    // support reference_wrappers around such types as well
    template <typename T>
    inline auto operator()(std::reference_wrapper<T> thing) const
        -> decltype(thing.get().geometryId(), Acts::GeometryIdentifier()) {
        return thing.get().geometryId();
    }
};

struct CompareGeometryId {
    // indicate that comparisons between keys and full objects are allowed.
    using is_transparent = void;
    // compare two elements using the automatic key extraction.
    template <typename Left, typename Right>
    constexpr bool operator()(Left&& lhs, Right&& rhs) const {
        return GeometryIdGetter()(lhs) < GeometryIdGetter()(rhs);
    }
};

/// A source link that stores just an index.
class IndexSourceLink final {
    public:
        using Index = uint32_t;

        /// Construct from geometry identifier and index.
        constexpr IndexSourceLink(Acts::GeometryIdentifier gid, Index idx)
            : m_geometryId(gid), m_index(idx) {}

        /// Construct an invalid source link. Must be default constructible to
        /// satisfy SourceLinkConcept.
        IndexSourceLink() = default;
        IndexSourceLink(const IndexSourceLink&) = default;
        IndexSourceLink(IndexSourceLink&&) = default;
        IndexSourceLink& operator=(const IndexSourceLink&) = default;
        IndexSourceLink& operator=(IndexSourceLink&&) = default;

        /// Access the index.
        constexpr Index index() const { return m_index; }
        
        Acts::GeometryIdentifier geometryId() const { return m_geometryId; }
        
        struct SurfaceAccessor {
            const Acts::Experimental::Detector& detector;
        
            const Acts::Surface* operator()(const Acts::SourceLink& sourceLink) const {
                const auto& indexSourceLink = sourceLink.get<IndexSourceLink>();
                return *detector.sensitiveHierarchyMap().find(
                    indexSourceLink.m_geometryId);
            }
        };

    private:
        Acts::GeometryIdentifier m_geometryId;
        Index m_index = 0;
        
        friend bool operator==(const IndexSourceLink& lhs,
            const IndexSourceLink& rhs) {
                return (lhs.geometryId() == rhs.geometryId()) &&
                    (lhs.m_index == rhs.m_index);
        }
        friend bool operator!=(const IndexSourceLink& lhs,
            const IndexSourceLink& rhs) {
                return !(lhs == rhs);
        }
};

/// Store elements that know their detector geometry id, e.g. simulation hits.
///
/// @tparam T type to be stored, must be compatible with `CompareGeometryId`
template <typename T>
using GeometryIdSet =
    boost::container::flat_set<T, CompareGeometryId>;

/// Container of index source links.
///
/// Since the source links provide a `.geometryId()` accessor, they can be
/// stored in an ordered geometry container.
using IndexSourceLinkContainer = GeometryIdSet<IndexSourceLink>;
