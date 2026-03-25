#pragma once

#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinUtility.hpp"

#include <unordered_map>

class NoMaterialDecorator : public Acts::IMaterialDecorator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// The surface binnings
    std::unordered_map<Acts::GeometryIdentifier, Acts::BinUtility>
        surfaceBinnings;
    /// Vetos for the material provider
    std::vector<Acts::GeometryIdentifier> vetos;
  };

  /// @brief Default destructor
  ~NoMaterialDecorator() override = default;

  /// @brief Constructor
  NoMaterialDecorator(const Config& cfg);

  /// @brief Decorate a surface
  ///
  /// @param surface the non-cost surface that is decorated
  void decorate(Acts::Surface& surface) const override;

  /// @brief Decorate a TrackingVolume
  ///
  /// @param volume the non-cost volume that is decorated
  void decorate(Acts::TrackingVolume& volume) const override;

 private:
  /// Private access to the configuration
  Config m_cfg;
};
