#pragma once

#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <memory>

#include "TrackingPipeline/Alignment/AlignmentStore.hpp"

/// @brief Class performing update of the global alignment parameters
class GlobalAlignmentTransformUpdater {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Alignment store instance
    std::shared_ptr<AlignmentStore> alignmentStore;
  };

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  /// @param level logging level
  GlobalAlignmentTransformUpdater(const Config& cfg,
                                  Acts::Logging::Level level);

  /// @brief perform update of the alignment parameters
  ///
  /// @param element detector element to update the transform of
  /// @param gctx current geometry context
  /// @param deltaTranslation x, y, z alignment delta
  /// @param deltaAngle theta_x, theta_y, theta_z alignment delta
  ///
  /// @return true if the alignment update succeded, false otherwise
  bool updateAlignmentParameters(Acts::DetectorElementBase* element,
                                 const Acts::GeometryContext& gctx,
                                 const Acts::Vector3& deltaTranslation,
                                 const Acts::Vector3& deltaAngles) const;

 protected:
  /// @brief access to logging instance
  const Acts::Logger& logger() const { return *m_logger; }

 private:
  /// Configuration
  Config m_cfg;

  /// Vector of active alignment dofs indices
  std::vector<std::size_t> m_activeIdxs;

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};
