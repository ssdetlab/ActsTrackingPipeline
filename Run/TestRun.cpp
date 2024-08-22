#include "ActsLUXEPipeline/LUXEROOTDataReader.hpp"
#include "ActsLUXEPipeline/E320Geometry.hpp"
#include "ActsLUXEPipeline/IdealSeeder.hpp"
#include "ActsLUXEPipeline/TrackFitter.hpp"
#include "ActsLUXEPipeline/ConstantBoundedField.hpp"
#include "ActsLUXEPipeline/Sequencer.hpp"
#include "ActsLUXEPipeline/ROOTFittedTrackWriter.hpp"
#include "ActsLUXEPipeline/E320MagneticField.hpp"
#include "ActsLUXEPipeline/QuadrupoleMagField.hpp"
#include "ActsLUXEPipeline/LxBFields.hpp"
#include "ActsLUXEPipeline/BinnedMagneticField.hpp"

#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"

#include <filesystem>
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
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
    E320Geometry::GeometryOptions gOpt;

    // --------------------------------------------------------------
    // Detector setup

    // Set the path to the gdml file
    // and the names of the volumes to be converted
    std::string gdmlPath =
            "/Users/alonlevi/CLionProjects/LUXEPipeline/Zstuff/ettgeom_magnet_pdc_tracker.gdml";
    std::vector<std::string> names{"OPPPSensitive"};

    // Build the detector
    auto trackerBP =
            E320Geometry::makeBlueprintE320(gdmlPath, names, gOpt);
    auto detector =
            E320Geometry::buildE320Detector(std::move(trackerBP), gctx, gOpt, {});

    Acts::ObjVisualization3D volumeObj;
    Acts::ViewConfig pConfig = Acts::s_viewSensitive;

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
    E320MagneticField dipoleFieldFunc(xParams,yParams,zParams,bxval);
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
    
    auto transformPos = [](const Acts::Vector3& pos) {
        E320Geometry::GeometryOptions dipoleOpt;
        for (int i=0;i<3;i++) {
            if (pos[i]<=-1*dipoleOpt.dipoleBounds[i] ||
                pos[i]>dipoleOpt.dipoleBounds[i]) {
                std::cout<<"Altered"<<std::endl;
                return Acts::Vector3{50,0,0};
            }
        }
        return pos;
    };
    auto dipolePos = [transformPos, gOpt](const Acts::Vector3& pos) {
        return transformPos(pos-gOpt.dipoleTranslation);
    };
    auto quad1Pos = [transformPos, gOpt](const Acts::Vector3& pos) {
        return transformPos(pos-gOpt.quad1Translation);
    };
    auto quad2Pos = [transformPos, gOpt](const Acts::Vector3& pos) {
        return transformPos(pos-gOpt.quad2Translation);
    };
    auto quad3Pos = [transformPos, gOpt](const Acts::Vector3& pos) {
        return transformPos(pos-gOpt.quad3Translation);
    };
// map (Bx,By,Bz) -> (Bx,By,Bz)
    auto transformBField = [](const Acts::Vector3& field, const Acts::Vector3&) {
        return field;
    };

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

    return 0;

}
