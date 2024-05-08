#pragma once

#include "ActsLUXEPipeline/LUXEGeometryConstraints.hpp"

#include "Acts/Detector/Blueprint.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

namespace LUXEGeometry {

/// @brief Make the blueprint for the LUXE detector
/// in the two arm configuration
///
/// @param gdmlPath path to the gdml file
/// @param names the names of the volumes to be converted
/// @param gOpt the geometry constraints
///
/// @return the Blueprint for the LUXE detector
std::unique_ptr<Acts::Experimental::Blueprint::Node> 
makeBlueprintLUXE(
    const std::string& gdmlPath,
    const std::vector<std::string>& names,
    const LUXEGeometry::GeometryOptions& gOpt);

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
        const LUXEGeometry::GeometryOptions& gOpt,
        const std::vector<Acts::GeometryIdentifier>& materialVetos);

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
        const LUXEGeometry::GeometryOptions& gOpt, 
        const std::string jsonMaterialPath,
        const std::vector<Acts::GeometryIdentifier>& materialVetos);


} // namespace LUXEGeometry
