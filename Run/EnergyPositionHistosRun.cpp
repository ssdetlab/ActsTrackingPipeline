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
#include <TFile.h>
#include <TTree.h>


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
            std::cout<<"Surface bounds: "<<surf->normal(gctx,surf->center(gctx),Acts::Vector3{0,1,0})<<std::endl;
        }
    }
    MeasurementResolutionMap resolutions = m;

    auto propagator = LUXENavigator::makePropagator<Acts::EigenStepper<>>(detector, BFieldPtr);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(2.65, 16.5);
// lower limit 0.32 GeV upper limit 3.4 GeV
    std::vector<LUXENavigator::Measurements> results;
    std::size_t sourceId = 1;
    for (int i=0;i<10;i++) {
        Acts::ActsScalar E = dis(gen);
        results.push_back(LUXENavigator::createMeasurements(propagator, gctx, magCtx,
                                                            LUXENavigator::makeParameters(E),
                                                            resolutions,sourceId));
        sourceId++;
        std::cout<<"Initial Energy : "<<E<<std::endl;
    };

    TFile *file = new TFile("HistogramData.root","RECREATE");
    TTree *tree = new TTree("tree", "TruthParameters");
    struct ROOTMeasurement {
        unsigned int id;
        unsigned int layer;
        float local_x;
        float local_y;
        float phi;
        float theta;
        float QOverP;
        float time;
    };

    ROOTMeasurement s;

    tree->Branch("id", &s.id, "id/I");
    tree->Branch("layer", &s.layer, "id/I");
    tree->Branch("local_x", &s.local_x);
    tree->Branch("local_y", &s.local_y);
    tree->Branch("phi", &s.phi);
    tree->Branch("theta", &s.theta);
    tree->Branch("QOverP", &s.QOverP);
    tree->Branch("time", &s.time);

    for (auto& result:results) {
        std::cout<<"globals size in results loop"<<result.globalPosition.size()<<std::endl;
        std::cout<<"geo IDs : "<<result.sourceLinks[0].m_geometryId<<std::endl;
        for (unsigned int l=0;l<result.truthParameters.size();l++) {
            s.id = result.eventId;
            s.layer = l+1;
            s.local_x = result.truthParameters[l][0];
            s.local_y = result.truthParameters[l][1];
            s.phi = result.truthParameters[l][2];
            s.theta = result.truthParameters[l][3];
            s.QOverP = result.truthParameters[l][4];
            s.time = result.truthParameters[l][5];
            tree->Fill();
//            if (l!=result.globalPosition.size()-1) {
//            Acts::GeometryView3D::drawArrowForward(
//                    volumeObj,result.globalPosition[l],
//                    result.globalPosition[l+1], 20, 20, pConfig);
//            }
        } //static_cast<float>
//        for (unsigned int l=1;l<result.fullTrack.size()-3;l++) {
//            if (result.fullTrack[l][1]<5000) {
//                Acts::GeometryView3D::drawSegment(
//                        volumeObj,result.fullTrack[l],
//                        result.fullTrack[l+1], pConfig);
//            }
//        }
    }

//    Acts::GeometryView3D::drawArrowForward(
//            volumeObj,Acts::Vector3{0,0,0},
//            Acts::Vector3{0,1200,0},1000,100, pConfig);
//    Acts::GeometryView3D::drawArrowForward(
//            volumeObj,Acts::Vector3{0,0,0},
//            Acts::Vector3{1200,0,0},1000,100, pConfig);
//    Acts::GeometryView3D::drawArrowForward(
//            volumeObj,Acts::Vector3{0,0,0},
//            Acts::Vector3{0,0,1200},1000,100, pConfig);

    file->Write();
    file->Close();
    delete file;
    volumeObj.write("volumes.obj");

    // Run all configured algorithms and return the appropriate status.
//    return sequencer.run();
      return 0;
} // main
