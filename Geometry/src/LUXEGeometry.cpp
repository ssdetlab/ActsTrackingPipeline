#include "ActsLUXEPipeline/LUXEGeometry.hpp"

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Plugins/Geant4/Geant4SurfaceProvider.hpp"

namespace LUXEGeometry {

std::tuple<G4VPhysicalVolume*, 
std::vector<std::string>, 
std::string>
MockupGeant4World() {
    const std::string gdmlPath = "two-arms-telescope.gdml";

    const int nLayers = 5;
    const int nChips = 5;
    const int nArms = 2;
    
    const double cellDimX = 0.4;
    const double cellDimY = 0.3;
    const double cellDimZ = 0.1;
    
    const double armOffset = 10;

    // Get nist material manager
    G4NistManager* nist = G4NistManager::Instance();

    // Option to switch on/off checking of volumes overlaps
    //
    G4bool checkOverlaps = true;
    
    //
    // World
    //
    G4double worldSizeXY = 50;
    G4double worldSizeZ = 50;
    G4Material* worldMat = nist->FindOrBuildMaterial("G4_Galactic");
    
    auto solidWorld = new G4Box("World", 0.5 * worldSizeXY, 0.5 * worldSizeXY,
                                0.5 * worldSizeZ);
    
    auto logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
    
    auto physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld,
                                        "World", nullptr, false, 0, checkOverlaps);
    
    G4Material* material = nist->FindOrBuildMaterial("G4_He");
    
    std::vector<std::string> names;
    for (int nArm = 0; nArm < nArms; nArm++) {
        for (int nLayer = 0; nLayer < nLayers; nLayer++) {
            for (int nChip = 0; nChip < nChips; nChip++) {
                int sign = (nArm == 0) ? 1 : -1;
                double posX = sign * (armOffset + nChip);
                double posY = 0;
                double posZ = nLayer;
                G4ThreeVector pos = G4ThreeVector(posX, posY, posZ);
                G4ThreeVector dims = G4ThreeVector(cellDimX, cellDimY, cellDimZ);
        
                int cellId = nChips * nLayers * nArm + nChips * nLayer + nChip;
                std::string name = "cell" + std::to_string(cellId);
        
                names.push_back(name);
        
                // Box cell
                auto solidCell = new G4Box(name, dims[0], dims[1], dims[2]);
        
                G4LogicalVolume* cellLogical =
                    new G4LogicalVolume(solidCell, material, "cellSensitive");
        
                new G4PVPlacement(nullptr, pos, cellLogical, name, logicWorld, false, 0,
                                true);
            }
        }
    }

    // Write GDML file
    G4GDMLParser parser;
    parser.SetOutputFileOverwrite(true);
    parser.Write(gdmlPath, physWorld);

    return std::make_tuple(physWorld, names, gdmlPath);
}

std::shared_ptr<const Acts::Experimental::Detector> 
    buildLUXEDetector(std::string gdmlPath, 
        std::vector<std::string> names, 
        Acts::GeometryContext gctx) {
            // Default template parameters are fine
            // when using names as identifiers
            auto spFullCfg = Acts::Experimental::Geant4SurfaceProvider<>::Config();
            spFullCfg.gdmlPath = gdmlPath;
            spFullCfg.surfacePreselector =
                std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(names,
                    true);
    
            auto spFull = std::make_shared<Acts::Experimental::Geant4SurfaceProvider<>>(
                spFullCfg, Acts::Experimental::Geant4SurfaceProvider<>::kdtOptions(),
                false);
    
            auto lbFullCfg = Acts::Experimental::LayerStructureBuilder::Config();
            lbFullCfg.surfacesProvider = spFull;
    
            auto lbFull =
                std::make_shared<Acts::Experimental::LayerStructureBuilder>(lbFullCfg);
    
            auto [sFull, vFull, suFull, vuFull] = lbFull->construct(gctx);
            std::cout<<sFull.size()<<std::endl;
//            assert(sFull.size() == names.size());

            // Create the detector
            auto bounds = std::make_unique<Acts::CuboidVolumeBounds>(3, 3, 3);
            auto volume = Acts::Experimental::DetectorVolumeFactory::construct(
                Acts::Experimental::defaultPortalAndSubPortalGenerator(), gctx,
                "volume", Acts::Transform3::Identity(), std::move(bounds), {}, {},
                Acts::Experimental::tryNoVolumes(),
                Acts::Experimental::tryAllPortalsAndSurfaces());
            volume->assignGeometryId(1);
            auto detector = Acts::Experimental::Detector::makeShared(
                "detector", {volume}, Acts::Experimental::tryRootVolumes());
            
            return detector;
}

} // namespace LUXEGeometry
