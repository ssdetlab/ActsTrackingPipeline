#pragma once

#include "Acts/Detector/Blueprint.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Plugins/Geant4/Geant4SurfaceProvider.hpp"
#include "ActsLUXEPipeline/LUXEGeometryConstraints.hpp"

#include <iostream>
#include <tuple>

namespace LUXEGeometry {

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
std::shared_ptr<Acts::Experimental::LayerStructureBuilder> 
makeLayerBuilder(
    const G4VPhysicalVolume* world,
    const G4Transform3D& transform,
    const std::vector<std::string>& names,
    const std::array<std::tuple<Acts::ActsScalar,
        Acts::ActsScalar>, kDim>& ranges,
    const std::array<Acts::BinningValue, kDim>& binningValues) {
        auto spCfg = 
            typename Acts::Experimental::Geant4SurfaceProvider<kDim>::Config();
        spCfg.g4World = world;
        spCfg.worldTransform = transform;
        spCfg.surfacePreselector =
            std::make_shared<
                Acts::Geant4PhysicalVolumeSelectors::NameSelector>(names, true);

        auto kdtDOpt = 
            typename Acts::Experimental::Geant4SurfaceProvider<kDim>::kdtOptions();
        for (std::size_t i = 0; i < kDim; i++) {
            kdtDOpt.range[i].set(
                std::get<0>(ranges[i]), std::get<1>(ranges[i]));
            kdtDOpt.binningValues[i] = binningValues[i];
        }

        auto sp = std::make_shared<
            Acts::Experimental::Geant4SurfaceProvider<kDim>>(spCfg, kdtDOpt);

        auto lbCfg = Acts::Experimental::LayerStructureBuilder::Config();
        lbCfg.surfacesProvider = sp;
        auto lb =
            std::make_shared<Acts::Experimental::LayerStructureBuilder>(lbCfg);
        return lb;
}

/// @brief Make the blueprint for the LUXE detector
/// in the positron arm configuration
///
/// @param gdmlPath path to the gdml file
/// @param names the names of the volumes to be converted
/// @param gOpt the geometry constraints
///
/// @return the Blueprint for the LUXE detector
std::unique_ptr<Acts::Experimental::Blueprint::Node> 
makeBlueprintPositron(
    const std::string& gdmlPath,
    const std::vector<std::string>& names,
    const GeometryOptions& gOpt);

/// @brief Build the LUXE detector
///
/// @param gdmlPath path to the gdml file
/// @param names the names of the volumes to be converted
/// @param gctx the geometry context
/// @param gOpt the geometry constraints
///
/// @return shared pointer to the detector object
std::shared_ptr<const Acts::Experimental::Detector> 
    buildLUXEDetector(
        const std::unique_ptr<Acts::Experimental::Blueprint::Node> 
            detectorBpr,
        const Acts::GeometryContext& gctx,
        const GeometryOptions& gOpt);

} // namespace LUXEGeometry
