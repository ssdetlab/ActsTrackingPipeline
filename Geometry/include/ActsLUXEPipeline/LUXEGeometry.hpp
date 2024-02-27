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

/// @brief Construct the test Geant4 world
/// of thwo arms telescope with 5 layers of 5 chips
///
/// @return a tuple of the world volume, 
/// the chips names and the gdml file path
std::tuple<G4VPhysicalVolume*, 
std::vector<std::string>, 
std::string>
MockupGeant4World();

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
