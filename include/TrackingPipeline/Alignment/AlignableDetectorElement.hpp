#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>

#include "TrackingPipeline/Alignment/AlignmentContext.hpp"

class AlignableDetectorElement : public Acts::DetectorElementBase {
 public:
  AlignableDetectorElement(const std::shared_ptr<Acts::Surface>& surface,
                           const Acts::Transform3& transform)
      : m_surface(surface), m_transform(transform) {}

  ~AlignableDetectorElement() = default;

  /// Return the transform for the Element proxy mechanism
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  const Acts::Transform3& transform(
      const Acts::GeometryContext& gctx) const override {
    if (!gctx.hasValue()) {
      return nominalTransform();
    }
    const auto& store = gctx.get<AlignmentContext&>().alignmentStore();
    const auto& geoId = m_surface->geometryId();
    if (store.contains(geoId)) {
      return store.at(geoId);
    } else {
      return nominalTransform();
    }
  }

  const Acts::Transform3& nominalTransform() const { return m_transform; }

  /// Return surface representation - const return pattern
  const Acts::Surface& surface() const override { return *m_surface; }

  /// Non-const return pattern
  Acts::Surface& surface() override { return *m_surface; }

  /// Returns the thickness of the module
  ///
  /// @return double that indicates the thickness of the module
  double thickness() const override {
    if (m_surface->surfaceMaterial() == nullptr) {
      return 0;
    }
    return m_surface->surfaceMaterial()
        ->materialSlab(m_surface->center(Acts::GeometryContext()))
        .thickness();
  }

 private:
  std::shared_ptr<Acts::Surface> m_surface;
  Acts::Transform3 m_transform;
};
