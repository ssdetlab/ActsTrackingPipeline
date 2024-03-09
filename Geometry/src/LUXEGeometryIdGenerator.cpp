// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsLUXEPipeline/LUXEGeometryIdGenerator.hpp"

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Surfaces/Surface.hpp"

Acts::Experimental::IGeometryIdGenerator::GeoIdCache
Acts::Experimental::LUXEGeometryIdGenerator::generateCache() const {
    return Cache{};
}

void Acts::Experimental::LUXEGeometryIdGenerator::assignGeometryId(
        IGeometryIdGenerator::GeoIdCache& cache, DetectorVolume& dVolume) const {
    auto& ccache = std::any_cast<Cache&>(cache);

    ACTS_VERBOSE("Processing volume " << dVolume.name());
    // Set to the volume itself
    if (dVolume.geometryId().volume() == 0 || m_cfg.overrideExistingIds) {
        ++ccache.volumeCount;
        GeometryIdentifier geoID = volumeId(ccache);
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

void Acts::Experimental::LUXEGeometryIdGenerator::assignGeometryId(
        IGeometryIdGenerator::GeoIdCache& cache, Portal& portal) const {
    auto& ccache = std::any_cast<Cache&>(cache);

    auto& pSurface = portal.surface();
    if (pSurface.geometryId().boundary() == 0 || m_cfg.overrideExistingIds) {
        GeometryIdentifier geoID = volumeId(ccache, false);
        geoID.setBoundary(++ccache.portalCount);
        ACTS_VERBOSE("Assigning portal id " << ccache.portalCount);
        pSurface.assignGeometryId(geoID);
    }
}

void Acts::Experimental::LUXEGeometryIdGenerator::assignGeometryId(
        IGeometryIdGenerator::GeoIdCache& cache, Surface& surface) const {
    auto& ccache = std::any_cast<Cache&>(cache);

    auto rGeoID = surface.geometryId();
    auto geoID = volumeId(ccache, false);
    if (!m_cfg.overrideExistingIds && rGeoID.value() != 0) {
        return;
    } else if ((rGeoID.sensitive() == 0 && rGeoID.passive() == 0) ||
               m_cfg.overrideExistingIds) {
        ACTS_VERBOSE("Assigning sensitive id " <<
            (ccache.volumeCount+std::floor(ccache.sensitiveCount/9)-2)*10+ccache.sensitiveCount%9+1);
        geoID.setSensitive(
            (ccache.volumeCount+std::floor(ccache.sensitiveCount/9)-2)*10+ccache.sensitiveCount%9+1);
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

Acts::GeometryIdentifier Acts::Experimental::LUXEGeometryIdGenerator::volumeId(
        Cache& cache, bool incrementLayer) const {
    GeometryIdentifier geoID(0u);
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