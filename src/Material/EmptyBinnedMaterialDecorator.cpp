#include "TrackingPipeline/Material/EmptyBinnedMaterialDecorator.hpp"

#include "Acts/Material/ProtoSurfaceMaterial.hpp"

EmptyBinnedMaterialDecorator::EmptyBinnedMaterialDecorator(const Config& cfg)
    : m_cfg(cfg) {}

void EmptyBinnedMaterialDecorator::decorate(Acts::Surface& surface) const {
  auto surfaceMaterial = std::make_shared<Acts::ProtoSurfaceMaterial>(
      m_cfg.surfaceBinnings.at(surface.geometryId()));
  surface.assignSurfaceMaterial(surfaceMaterial);
};

void EmptyBinnedMaterialDecorator::decorate(
    Acts::TrackingVolume& volume) const {};
