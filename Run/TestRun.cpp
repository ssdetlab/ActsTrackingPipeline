#include "ActsLUXEPipeline/LUXEROOTDataReader.hpp"
#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/IdealSeeder.hpp"
#include "ActsLUXEPipeline/TrackFitter.hpp"
#include "ActsLUXEPipeline/ConstantBoundedField.hpp"
#include "ActsLUXEPipeline/Sequencer.hpp"
#include "ActsLUXEPipeline/ROOTFittedTrackWriter.hpp"
#include "ActsLUXEPipeline/LxBFields.hpp"

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
using namespace CLHEP;

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
        "/Users/alonlevi/CLionProjects/LUXEPipeline/Zstuff/ettgeom_magnet_pdc_tracker.gdml";
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

    G4double bxval, byval, bzval;
    bxval = byval = bzval = 0.0;

    byval = -0.95*tesla;

    std::map<std::string, std::vector<std::string>> model_params {
            {{"x"}, {"f_fd", "-165.0 165.0 7.7 7.7 mm"}},
            {{"y"}, {"const", "-30.0 30.0 mm"}},
            {{"z"}, {"f_fd", "-616.0 622.0 28.66 28.91 mm"}}
    };

    std::vector<FieldDistribution*> bcj(9,0);
    bcj[3] = LxBFieldAux::CreateFieldDistribution(model_params["x"][0], model_params["x"][1]);
    bcj[4] = LxBFieldAux::CreateFieldDistribution(model_params["y"][0], model_params["y"][1]);
    bcj[5] = LxBFieldAux::CreateFieldDistribution(model_params["z"][0], model_params["z"][1]);
    std::for_each(bcj.begin(), bcj.end(), [](FieldDistribution* &bd) {if (!bd) bd = new FieldDistribution(); });

    typedef VFieldComponentDistrib<FieldDistribution, FieldDistribution, FieldDistribution> LxTFieldComponent;
    LxTFieldComponent *bx = new LxTFieldComponent(bcj[0], bcj[1], bcj[2], bxval);
    LxTFieldComponent *by = new LxTFieldComponent(bcj[3], bcj[4], bcj[5], byval);
    LxTFieldComponent *bz = new LxTFieldComponent(bcj[6], bcj[7], bcj[8], bzval);

    double IPMagnetXpos = 0.0 *mm;
    double IPMagnetYpos = 0.0 *mm;
    double IPMagnetZpos = (150.0 + 55.0) *cm;

    G4ThreeVector magpos(IPMagnetXpos, IPMagnetYpos, IPMagnetZpos);
    LxBField *bfield = new LxBField(bx, by, bz, magpos);

//////////////////////////////////////////////////////////////////
    double MagField[6];
    double ppos[4] = {0.0, 0.0, IPMagnetZpos, 0.0};
    bfield->GetFieldValue(ppos, MagField);

    volumeObj.write("volumes.obj");

    return 0;

}
