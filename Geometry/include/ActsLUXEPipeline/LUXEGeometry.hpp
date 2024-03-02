#pragma once

#include "Acts/Detector/Detector.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

#include "G4VPhysicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GDMLParser.hh"

#include <iostream>

namespace LUXEGeometry {

/// @brief Build the LUXE detector
///
/// @param gdmlPath path to the gdml file
/// @param names the names of the volumes to be converted
///
/// @return shared pointer to the detector object
std::shared_ptr<const Acts::Experimental::Detector> 
    buildLUXEDetector(std::string gdmlPath, 
        std::vector<std::string> names, 
        Acts::GeometryContext gctx);

} // namespace LUXEGeometry
