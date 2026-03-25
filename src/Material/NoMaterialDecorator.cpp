#include "TrackingPipeline/Material/NoMaterialDecorator.hpp"

#include "Acts/Material/ProtoSurfaceMaterial.hpp"

NoMaterialDecorator::NoMaterialDecorator(const Config& cfg) : m_cfg(cfg) {}

void NoMaterialDecorator::decorate(Acts::Surface& surface) const {
  if (std::find(m_cfg.vetos.begin(), m_cfg.vetos.end(), surface.geometryId()) !=
      m_cfg.vetos.end()) {
    return;
  }
  auto surfaceMaterial = std::make_shared<Acts::ProtoSurfaceMaterial>(
      m_cfg.surfaceBinnings.at(surface.geometryId()));
  surface.assignSurfaceMaterial(surfaceMaterial);
};

void NoMaterialDecorator::decorate(Acts::TrackingVolume& volume) const {};
