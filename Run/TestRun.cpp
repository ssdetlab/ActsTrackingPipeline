#include "ActsLUXEPipeline/LUXEROOTDataReader.hpp"
#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/IdealSeeder.hpp"
#include "ActsLUXEPipeline/TrackFitter.hpp"
#include "ActsLUXEPipeline/ConstantBoundedField.hpp"
#include "ActsLUXEPipeline/Sequencer.hpp"
#include "ActsLUXEPipeline/ROOTFittedTrackWriter.hpp"
<<<<<<< Updated upstream
=======
#include "ActsLUXEPipeline/DipoleMagField.hpp"
#include "ActsLUXEPipeline/QuadrupoleMagField.hpp"
#include "ActsLUXEPipeline/LxBFields.hpp"
#include "ActsLUXEPipeline/BinnedMagneticField.hpp"
>>>>>>> Stashed changes

#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"

#include <filesystem>

#include "Acts/Plugins/Geant4/Geant4SurfaceProvider.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"

#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

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

<<<<<<< Updated upstream
    auto spFull = std::make_shared<Acts::Experimental::Geant4SurfaceProvider<>>(
        spFullCfg, Acts::Experimental::Geant4SurfaceProvider<>::kdtOptions());
=======
    for (auto& v : detector->volumes()) {
        std::cout << v->name() << std::endl;
        if (v->name()=="Dipole" || v->name()=="Tracker_gap_0") {
            std::cout << v->volumeBounds().values()[0] << std::endl;
            std::cout << v->center(gctx) << std::endl;
        }
        Acts::GeometryView3D::drawDetectorVolume(volumeObj,
                                                 *v,gctx);
        for (auto& s : v->surfaces()) {
            std::cout << s->center(gctx).transpose() << "   " << s->normal(gctx, s->center(gctx), Acts::Vector3(1,0,0)).transpose() << std::endl;
        }
    }


    Acts::ActsScalar bxval = 0.31_T;
    std::vector<double> GridLimits =
            detector->findDetectorVolume("Dipole")->volumeBounds().values();

    Acts::Vector2 xParams(-30.0,30.0);
    Acts::Vector4 yParams(-165.0,165.0,7.7,7.7);
    Acts::Vector4 zParams(-457.0,457.0,25.0,25.0);
    DipoleMagField dipoleFieldFunc(xParams,yParams,zParams,bxval);
    Acts::ActsScalar quad1grad = 4_T / (1_m);
    Acts::ActsScalar quad2grad = -7_T / (1_m);
    Acts::ActsScalar quad3grad = 4_T / (1_m);
    QuadrupoleMagField quad1FieldFunc(quad1grad);
    QuadrupoleMagField quad2FieldFunc(quad2grad);
    QuadrupoleMagField quad3FieldFunc(quad3grad);

    Acts::MagneticFieldProvider::Cache dipoleCache = dipoleFieldFunc.makeCache(mctx);
    Acts::MagneticFieldProvider::Cache quad1Cache = quad1FieldFunc.makeCache(mctx);
    Acts::MagneticFieldProvider::Cache quad2Cache = quad2FieldFunc.makeCache(mctx);
    Acts::MagneticFieldProvider::Cache quad3Cache = quad3FieldFunc.makeCache(mctx);
>>>>>>> Stashed changes
    
    auto lbFullCfg = Acts::Experimental::LayerStructureBuilder::Config();
    lbFullCfg.surfacesProvider = spFull;
    
    auto lbFull =
        std::make_shared<Acts::Experimental::LayerStructureBuilder>(lbFullCfg);
    
    auto [sFull, vFull, suFull, vuFull] = lbFull->construct(gctx);

<<<<<<< Updated upstream
    int k = 0;
    for (auto& surf : sFull) {
        std::cout << k << " Surface: " << surf->center(gctx).transpose() << std::endl;
        k++;
    }
=======
    vGridOptions gridOpt;
    gridOpt.xBins = {-1*GridLimits[0],-200,-166, 0,0.5, 164,165,166,GridLimits[0]};
    gridOpt.yBins = {-1*GridLimits[1],-1, 0,0.5, 28,29,30,51,52,GridLimits[1]};
    gridOpt.zBins = {-1*GridLimits[2],-2,-1, 0,0.5, 1.0,2,GridLimits[2]};

    auto DipoleField = buildBinnedBField(dipoleFieldFunc, dipolePos, transformBField, gridOpt, mctx);
    auto DipoleFieldPtr = std::make_shared<Acts::InterpolatedBFieldMap<vGrid>>(DipoleField);
    auto Quad1Field = buildBinnedBField(quad1FieldFunc, quad1Pos, transformBField, gridOpt, mctx);
    auto Quad1FieldPtr = std::make_shared<Acts::InterpolatedBFieldMap<vGrid>>(Quad1Field);
    auto Quad2Field = buildBinnedBField(quad2FieldFunc, quad2Pos, transformBField, gridOpt, mctx);
    auto Quad2FieldPtr = std::make_shared<Acts::InterpolatedBFieldMap<vGrid>>(Quad2Field);
    auto Quad3Field = buildBinnedBField(quad3FieldFunc, quad3Pos, transformBField, gridOpt, mctx);
    auto Quad3FieldPtr = std::make_shared<Acts::InterpolatedBFieldMap<vGrid>>(Quad3Field);
    
    
//    std::cout<<"FieldPtr"<<(DipoleFieldPtr->getField(Acts::Vector3{0,0,13140})).value()<<std::endl;
//    std::cout<<bxval<<std::endl;
    volumeObj.write("E320volumes.obj");
>>>>>>> Stashed changes

    return 0;

}
