#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>
#include <stdexcept>

#include "TrackingPipeline/Alignment/AlignmentContext.hpp"

class AlignableDetectorElement : public Acts::DetectorElementBase {
 public:
  AlignableDetectorElement(std::shared_ptr<Acts::Surface> surface,
                           const Acts::Transform3& transform)
      : m_surface(surface), m_transform(transform) {}

  ~AlignableDetectorElement() = default;

  /// Return the transform for the Element proxy mechanism
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  const Acts::Transform3& transform(
      const Acts::GeometryContext& gctx) const override {
    if (!gctx.hasValue()) {
      return m_transform;
    }
    auto alignContext = gctx.get<AlignmentContext>();
    if (alignContext.alignmentStore != nullptr &&
        alignContext.alignmentIndex < 2) {
      return (*(alignContext.alignmentStore))[alignContext.alignmentIndex];
    } else {
      throw std::runtime_error("Invalid alignment context");
    }
  }

  /// Return surface representation - const return pattern
  const Acts::Surface& surface() const override { return *m_surface; }

  /// Non-const return pattern
  Acts::Surface& surface() override { return *m_surface; }

  /// Returns the thickness of the module
  /// @return double that indicates the thickness of the module
  ///
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
