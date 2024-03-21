#include "ActsLUXEPipeline/Sequencer.hpp"

#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXEGeometryConstraints.hpp"
#include "ActsLUXEPipeline/LUXEMagneticField.hpp"
#include "ActsLUXEPipeline/LUXENavigator.hpp"
#include "ActsLUXEPipeline/LUXEMeasurementsCreator.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"

#include <filesystem>
#include <string>
#include <iostream>
#include <random>

#include "ActsLUXEPipeline/Utils.hpp"

/// @brief Run the propagation through
/// a uniform energy spectrum and record the
/// energy vs position histograms for each layer
int main() {
    using namespace LUXENavigator;
    Acts::Logging::Level logLevel = Acts::Logging::VERBOSE;

    // setup the sequencer first w/ config derived from options
    Sequencer::Config seqCfg;
    seqCfg.events = 10;
    seqCfg.numThreads = -1;
    Sequencer sequencer(seqCfg);

//    LUXEROOTReader::LUXEROOTSimDataReader::Config readerCfg
//        = LUXEROOTReader::defaultSimConfig();
//    readerCfg.dataCollection = "SourceLink";
//    std::string pathToDir = "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_dataInRootFormat/SignalNextTrial_e1gpc_10.0_1";
    // map (x,y,z) -> (x,y,z)

    auto transformPos = [](const Acts::Vector3& pos) {
        LUXEGeometry::GeometryOptions gOpt;
        for (int i=0;i<3;i++) {
            if (pos[i]<gOpt.MagneticFieldBounds[i].first ||
                pos[i]>gOpt.MagneticFieldBounds[i].second) {
                return Acts::Vector3{0,-1,0};
            }
        }
        return pos;
    };

    // map (Bx,By,Bz) -> (Bx,By,Bz)
    auto transformBField = [](const Acts::Vector3& field, const Acts::Vector3&) {
        return field;
    };

    LUXEMagneticField::GridOptions gridOpt;
    gridOpt.bins = {14u, 1400u, 2u};
    gridOpt.limits = {std::make_pair(-1,1300),
                      std::make_pair(-1,1300),
                      std::make_pair(-1,1300)};

    auto BField = LUXEMagneticField::buildLUXEBField(transformPos, transformBField, gridOpt);
    auto BFieldPtr = std::make_shared<LUXEMagneticField::BField_t>(BField);

    // Build the LUXE detector
    std::string gdmlPath = "lxgeomdump_ip_tracker_positron.gdml";
    std::vector<std::string> staves = {"OPPPSensitive"};
    std::vector<std::string> chamber = {"IPMagnetField","VCWindowPanel"};
    Acts::GeometryContext gctx;
    Acts::MagneticFieldContext magCtx;
    LUXEGeometry::GeometryOptions gOpt;
    auto positronArmBpr = LUXEGeometry::makeBlueprintPositron(gdmlPath, staves, gOpt);
    auto magneticChamberBpr = LUXEGeometry::makeBlueprintMagneticChamber(gdmlPath, chamber, gOpt);
    positronArmBpr->add(std::move(magneticChamberBpr));
    auto detector = LUXEGeometry::buildLUXEDetector(std::move(positronArmBpr), gctx, gOpt);

    MeasurementResolution resPixel = {MeasurementType::eLoc01,
                                      {gOpt.chipSizeX,
                                       gOpt.chipSizeY}};
    std::vector<std::pair<Acts::GeometryIdentifier,MeasurementResolution>> m;
    Acts::ViewConfig pConfig = Acts::s_viewSensitive;
    Acts::ObjVisualization3D volumeObj;
    for (auto& vol : detector->rootVolumes()) {
        std::cout<<"Surfaces size: "<<vol->surfaces().size()<<std::endl;
        std::cout<<"Volume Bounds: "<<vol->volumeBounds()<<std::endl;
        std::cout<<"Volume Transformation: "<<vol->transform().translation()<<std::endl;
//        Acts::GeometryView3D::drawDetectorVolume(
//                volumeObj, *(vol), gctx,
//                Acts::Transform3::Identity(), pConfig);
        for (auto& surf : vol->surfaces()) {
            std::cout<<"Assigning resolution to surface ID: "<<surf->geometryId()<<std::endl;
            if (vol->geometryId().volume()!=1) {
                Acts::GeometryView3D::drawSurface(
                        volumeObj, *(surf), gctx,
                        Acts::Transform3::Identity(), pConfig);
            }
            m.push_back(std::make_pair(surf->geometryId(),resPixel));
            std::cout<<"Surface x transform: "<<surf->center(gctx)[0]<<std::endl;
            std::cout<<"Surface y transform: "<<surf->center(gctx)[1]<<std::endl;
            std::cout<<"Surface z transform: "<<surf->center(gctx)[2]<<std::endl;
            std::cout<<"Surface bounds: "<<surf->bounds()<<std::endl;
        }
    }
    MeasurementResolutionMap resolutions = m;

    auto propagator = LUXENavigator::makePropagator<Acts::EigenStepper<>>(detector, BFieldPtr);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.5, 16.5);
// lower limit 0.32 GeV upper limit 3.4 GeV

//    std::vector<Acts::ActsScalar> test_E{0.1,0.3,0.4,0.5,0.6};

    std::vector<LUXENavigator::Measurements> results;
    std::size_t sourceId = 1;
    for (int i=0;i<100000;i++) {
        Acts::ActsScalar E = dis(gen);
        std::cout<<"Initial Energy : "<<E<<std::endl;
        results.push_back(LUXENavigator::createMeasurements(propagator, gctx, magCtx,
                                                            LUXENavigator::makeParameters(E),
                                                            resolutions,sourceId));
        sourceId++;
    };

    for (auto& result:results) {
        for (unsigned int l=0;l<result.truthParameters.size();l++) {
            if (l!=result.globalPosition.size()-1) {
            Acts::GeometryView3D::drawArrowForward(
                    volumeObj,result.globalPosition[l],
                    result.globalPosition[l+1], 20, 20, pConfig);
            }
        } //static_cast<float>
//        for (unsigned int l=1;l<result.fullTrack.size()-3;l++) {
//            if (result.fullTrack[l][1]<6000) {
//                Acts::GeometryView3D::drawSegment(
//                        volumeObj,result.fullTrack[l],
//                        result.fullTrack[l+1], pConfig);
//            }
//        }
    }

    std::string filename = "hist_data.root";
    HistogramDatawriter(results,filename);
    volumeObj.write("volumes.obj");



    // Run all configured algorithms and return the appropriate status.
//    return sequencer.run();
      return 0;
} // main
