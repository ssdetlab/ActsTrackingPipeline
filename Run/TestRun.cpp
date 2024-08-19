#include "ActsLUXEPipeline/LUXEROOTDataReader.hpp"
#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/IdealSeeder.hpp"
#include "ActsLUXEPipeline/TrackFitter.hpp"
#include "ActsLUXEPipeline/ConstantBoundedField.hpp"
#include "ActsLUXEPipeline/Sequencer.hpp"
#include "ActsLUXEPipeline/ROOTFittedTrackWriter.hpp"

#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"

#include <filesystem>

#include "Acts/Plugins/Geant4/Geant4SurfaceProvider.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"

#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"

using ActionList = Acts::ActionList<>;
using AbortList = Acts::AbortList<Acts::EndOfWorldReached>;

using Propagator = Acts::Propagator<
    Acts::EigenStepper<>, 
    Acts::Experimental::DetectorNavigator>;
using PropagatorOptions =
    typename Propagator::template Options<ActionList, AbortList>;

using Trajectory = Acts::VectorMultiTrajectory;
using TrackContainer = Acts::VectorTrackContainer;
using KF = Acts::KalmanFitter<Propagator, Trajectory>;

using namespace Acts::UnitLiterals;

int main() {
    // Set the log level
    Acts::Logging::Level logLevel = Acts::Logging::VERBOSE;

    // Dummy context and options
    Acts::GeometryContext gctx;
    Acts::MagneticFieldContext mctx;
    Acts::CalibrationContext cctx;
    LUXEGeometry::GeometryOptions gOpt;

    // --------------------------------------------------------------
    // LUXE detector setup

    // Set the path to the gdml file
    // and the names of the volumes to be converted
    std::string gdmlPath = 
        "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/E320Pipeline_gdmls/ettgeom_magnet_pdc_tracker.gdml";
    std::vector<std::string> names{"OPPPSensitive"};

    /// Read the gdml file and get the world volume
    G4GDMLParser parser;
    parser.Read(gdmlPath, false);
    auto world = parser.GetWorldVolume();

    // Default template parameters are fine
    // when using names as identifiers
    auto spFullCfg = Acts::Experimental::Geant4SurfaceProvider<>::Config();
    spFullCfg.g4World = world;
    spFullCfg.surfacePreselector =
        std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(names,
                                                                            true);

    auto spFull = std::make_shared<Acts::Experimental::Geant4SurfaceProvider<>>(
        spFullCfg, Acts::Experimental::Geant4SurfaceProvider<>::kdtOptions());
    
    auto lbFullCfg = Acts::Experimental::LayerStructureBuilder::Config();
    lbFullCfg.surfacesProvider = spFull;
    
    auto lbFull =
        std::make_shared<Acts::Experimental::LayerStructureBuilder>(lbFullCfg);
    
    auto [sFull, vFull, suFull, vuFull] = lbFull->construct(gctx);

    Acts::ObjVisualization3D volumeObj;
    Acts::ViewConfig pConfig = Acts::s_viewSensitive;

    int k = 0;
    for (auto& surf : sFull) {
        std::cout << k << " Surface: " << surf->center(gctx).transpose() << std::endl;
        k++;
        Acts::GeometryView3D::drawSurface(
                volumeObj,*(surf),gctx,
                Acts::Transform3::Identity(),pConfig);
    }

    volumeObj.write("volumes.obj");

    return 0;

}
