#pragma once

#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinUtility.hpp"

#include <unordered_map>

/// @brief Material decorator assigning surfaces with binned empty material
///
/// The decorator assignes binned material to the surfaces with empty material 
/// contents. This binned material is further filled by the material mapping 
/// algorithm. The binned material decoration is requred since the default 
/// geometry construction assumes homogenous surface material
class EmptyBinnedMaterialDecorator : public Acts::IMaterialDecorator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// The surface binnings
    std::unordered_map<Acts::GeometryIdentifier, Acts::BinUtility>
        surfaceBinnings;
  };

  /// @brief Default destructor
  ~EmptyBinnedMaterialDecorator() override = default;

  /// @brief Constructor
  explicit EmptyBinnedMaterialDecorator(const Config& cfg);

  /// @brief Decorate a surface
  ///
  /// @param surface the non-cost surface that is decorated
  void decorate(Acts::Surface& surface) const override;

  /// @brief Decorate a TrackingVolume
  ///
  /// @param volume the non-cost volume that is decorated
  void decorate(Acts::TrackingVolume& volume) const override;

 private:
  /// Configuration
  Config m_cfg;
};
