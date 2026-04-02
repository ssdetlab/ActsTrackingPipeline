#pragma once

#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <memory>

#include "TrackingPipeline/Alignment/AlignmentStore.hpp"

class GlobalAlignmentTransformUpdater {
 public:
  struct Config {
    std::shared_ptr<AlignmentStore> alignmentStore;
  };

  GlobalAlignmentTransformUpdater(const Config& cfg,
                                  Acts::Logging::Level level);

  bool updateAlignmentParameters(Acts::DetectorElementBase* element,
                                 const Acts::GeometryContext& gctx,
                                 const Acts::Vector3& deltaTranslation,
                                 const Acts::Vector3& deltaAngles) const;

 protected:
  const Acts::Logger& logger() const { return *m_logger; }

 private:
  Config m_cfg;

  std::vector<std::size_t> m_activeIdxs;

  std::unique_ptr<const Acts::Logger> m_logger;
};
