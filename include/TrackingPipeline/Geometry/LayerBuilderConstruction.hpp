#pragma once

#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Plugins/Geant4/Geant4SurfaceProvider.hpp"

#include <tuple>

/// @brief Create internal structure builder
/// for the detector
///
/// @tparam kDim the dimension for the detector binning
///
/// @param world the G4 world volume
/// @param transform the transformation to apply to the world
/// @param names the names of the G4 volumes to be converted
/// @param ranges the ranges for the binning
/// @param binningValues axes to be binned
///
/// @return InternalStructureBuilder for the selected surfaces
template <std::size_t kDim = 1u>
std::shared_ptr<Acts::Experimental::LayerStructureBuilder> makeLayerBuilder(
    const G4VPhysicalVolume* world, const G4Transform3D& transform,
    const std::vector<std::string>& names,
    const std::array<std::tuple<double, double>, kDim>& ranges,
    const std::array<Acts::BinningValue, kDim>& binningValues,
    bool convertMaterial = false) {
  auto spCfg =
      typename Acts::Experimental::Geant4SurfaceProvider<kDim>::Config();
  spCfg.g4World = world;
  spCfg.worldTransform = transform;
  spCfg.surfacePreselector =
      std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(names,
                                                                          true);
  spCfg.convertMaterial = convertMaterial;

  auto kdtDOpt =
      typename Acts::Experimental::Geant4SurfaceProvider<kDim>::kdtOptions();
  for (std::size_t i = 0; i < kDim; i++) {
    kdtDOpt.range[i].set(std::get<0>(ranges[i]), std::get<1>(ranges[i]));
    kdtDOpt.binningValues[i] = binningValues[i];
  }

  auto sp = std::make_shared<Acts::Experimental::Geant4SurfaceProvider<kDim>>(
      spCfg, kdtDOpt);

  auto lbCfg = Acts::Experimental::LayerStructureBuilder::Config();
  lbCfg.surfacesProvider = sp;

  auto lb = std::make_shared<Acts::Experimental::LayerStructureBuilder>(lbCfg);
  return lb;
}