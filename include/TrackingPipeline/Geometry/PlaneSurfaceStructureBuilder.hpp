#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/RectangleBounds.hpp>

#include <memory>
#include <string>

class PlaneSurfaceStructureBuilder
    : public Acts::Experimental::IInternalStructureBuilder {
 public:
  struct Config {
    Acts::Transform3 transform;
    double halfX;
    double halfY;
  };

  PlaneSurfaceStructureBuilder(
      const Config& config,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "PlaneSurfaceStructureBuilder", Acts::Logging::INFO))
      : m_cfg(config) {};

  Acts::Experimental::InternalStructure construct(
      const Acts::GeometryContext& gctx) const final {
    Acts::Experimental::InternalStructure structure;
    structure.surfaces = {Acts::Surface::makeShared<Acts::PlaneSurface>(
        m_cfg.transform,
        std::make_shared<Acts::RectangleBounds>(m_cfg.halfX, m_cfg.halfY))};
    structure.surfacesUpdater = Acts::Experimental::tryAllPortalsAndSurfaces();
    structure.volumeUpdater = Acts::Experimental::tryNoVolumes();
    return structure;
  };

 private:
  Config m_cfg;
  /// Private access method to the logger
  const Acts::Logger& logger() const { return *m_logger; }

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};
