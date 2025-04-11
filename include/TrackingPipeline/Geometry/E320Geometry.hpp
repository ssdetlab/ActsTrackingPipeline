#pragma once

#include "Acts/Detector/Blueprint.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"

namespace E320Geometry {

/// @brief Make the blueprint for the E320 detector
///
/// @param gdmlPath path to the gdml file
/// @param names the names of the volumes to be converted
/// @param gOpt the geometry constraints
///
/// @return the Blueprint for the E320 detector
std::unique_ptr<Acts::Experimental::Blueprint::Node> makeBlueprintE320(
    const std::string& gdmlPath, const std::vector<std::string>& names,
    const E320Geometry::GeometryOptions& gOpt);

/// @brief Build the E320 detector
///
/// @param gdmlPath path to the gdml file
/// @param names the names of the volumes to be converted
/// @param gctx the geometry context
/// @param gOpt the geometry constraints
///
/// @return shared pointer to the detector object
std::shared_ptr<const Acts::Experimental::Detector> buildE320Detector(
    const std::unique_ptr<Acts::Experimental::Blueprint::Node> detectorBpr,
    const Acts::GeometryContext& gctx,
    const E320Geometry::GeometryOptions& gOpt,
    const std::vector<Acts::GeometryIdentifier>& materialVetos);

/// @brief Build the E320 detector
///
/// @param gdmlPath path to the gdml file
/// @param names the names of the volumes to be converted
/// @param gctx the geometry context
/// @param gOpt the geometry constraints
///
/// @return shared pointer to the detector object
std::shared_ptr<const Acts::Experimental::Detector> buildE320Detector(
    const std::unique_ptr<Acts::Experimental::Blueprint::Node> detectorBpr,
    const Acts::GeometryContext& gctx,
    const E320Geometry::GeometryOptions& gOpt,
    const std::string jsonMaterialPath,
    const std::vector<Acts::GeometryIdentifier>& materialVetos);

}  // namespace E320Geometry
