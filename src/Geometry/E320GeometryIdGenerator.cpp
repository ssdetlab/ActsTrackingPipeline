#include "TrackingPipeline/Geometry/E320GeometryIdGenerator.hpp"

namespace E320Geometry {

Acts::Experimental::IGeometryIdGenerator::GeoIdCache
E320GeometryIdGenerator::generateCache() const {
  return Cache{};
}

void E320GeometryIdGenerator::assignGeometryId(
    Acts::Experimental::IGeometryIdGenerator::GeoIdCache& cache,
    Acts::Experimental::DetectorVolume& dVolume) const {
  auto& ccache = std::any_cast<Cache&>(cache);

  ACTS_VERBOSE("Processing volume " << dVolume.name());
  // Set to the volume itself
  if (dVolume.geometryId().volume() == 0 || m_cfg.overrideExistingIds) {
    ++ccache.volumeCount;
    Acts::GeometryIdentifier geoID = volumeId(ccache);
    ACTS_VERBOSE("Assigning volume id " << geoID.volume());
    dVolume.assignGeometryId(geoID);
  }

  // Portals
  std::for_each(dVolume.portalPtrs().begin(), dVolume.portalPtrs().end(),
                [&](auto& portal) { assignGeometryId(cache, *portal); });

  // Surfaces
  std::for_each(dVolume.surfacePtrs().begin(), dVolume.surfacePtrs().end(),
                [&](auto& surface) { assignGeometryId(cache, *surface); });

  if (m_cfg.resetSubCounters) {
    ccache.portalCount = 0u;
    ccache.sensitiveCount = 0u;
    ccache.passiveCount = 0u;
  }

  // Sub volumes
  std::for_each(dVolume.volumePtrs().begin(), dVolume.volumePtrs().end(),
                [&](auto& volume) { assignGeometryId(cache, *volume); });
}

void E320GeometryIdGenerator::assignGeometryId(
    Acts::Experimental::IGeometryIdGenerator::GeoIdCache& cache,
    Acts::Experimental::Portal& portal) const {
  auto& ccache = std::any_cast<Cache&>(cache);

  auto& pSurface = portal.surface();
  if (pSurface.geometryId().boundary() == 0 || m_cfg.overrideExistingIds) {
    Acts::GeometryIdentifier geoID = volumeId(ccache, false);
    geoID.setBoundary(++ccache.portalCount);
    ACTS_VERBOSE("Assigning portal id " << ccache.portalCount);
    pSurface.assignGeometryId(geoID);
  }
}

void E320GeometryIdGenerator::assignGeometryId(
    Acts::Experimental::IGeometryIdGenerator::GeoIdCache& cache,
    Acts::Surface& surface) const {
  auto& ccache = std::any_cast<Cache&>(cache);

  auto rGeoID = surface.geometryId();
  auto geoID = Acts::GeometryIdentifier(0u);
  if (!m_cfg.overrideExistingIds && rGeoID.value() != 0) {
    return;
  } else if ((rGeoID.sensitive() == 0 && rGeoID.passive() == 0) ||
             m_cfg.overrideExistingIds) {
    Acts::Vector3 center = surface.center(Acts::GeometryContext());
    ACTS_VERBOSE("Processing surface " << center.transpose());
    std::int32_t geoIDval = 0u;

    // Determine the layer
    std::int32_t chipId = -1;
    for (auto [id, z] : m_cfg.gOpt.staveZ) {
      // These are already rotated surfaces
      if (std::abs(center.y() - z) < 1e-1) {
        chipId = id;
        break;
      }
    }
    if (chipId == -1) {
      geoID.setPassive(++ccache.passiveCount);
      ACTS_VERBOSE("No layer. Assigning passive id " << geoID.passive());
      surface.assignGeometryId(geoID);
      return;
    }
    geoIDval = chipId + 1;

    ACTS_VERBOSE("Assigning sensitive id " << geoIDval);
    geoID.setSensitive(geoIDval);
    ccache.sensitiveCount++;

    surface.assignGeometryId(geoID);
  } else if (rGeoID.sensitive() != 0 || rGeoID.passive() != 0) {
    ACTS_VERBOSE(
        "Surface already has a geometry id, only setting volume and layer id.");
    rGeoID.setVolume(geoID.volume());
    rGeoID.setLayer(geoID.layer());
    surface.assignGeometryId(rGeoID);
  }
}

Acts::GeometryIdentifier E320GeometryIdGenerator::volumeId(
    Cache& cache, bool incrementLayer) const {
  Acts::GeometryIdentifier geoID(0u);
  if (!m_cfg.containerMode) {
    geoID.setVolume(cache.volumeCount);
  } else {
    geoID.setVolume(m_cfg.containerId);
    if (incrementLayer) {
      ++cache.layerCount;
    }
    geoID.setLayer(cache.layerCount);
    ACTS_VERBOSE("Container mode: assigning volume id "
                 << m_cfg.containerId << ", layer id " << cache.layerCount);
  }
  return geoID;
}

}  // namespace E320Geometry
