#include "ActsLUXEPipeline/LUXEROOTDataReader.hpp"
#include "ActsLUXEPipeline/LUXESeeder.hpp"
#include "ActsLUXEPipeline/Sequencer.hpp"

#include <filesystem>

#include "ActsLUXEPipeline/LUXEGeometry.hpp"
#include "ActsLUXEPipeline/LUXEGeometryConstraints.hpp"
#include "ActsLUXEPipeline/LUXEMagneticField.hpp"
#include "ActsLUXEPipeline/LUXENavigator.hpp"
#include "ActsLUXEPipeline/LUXEMeasurementsCreator.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include <string>
#include <iostream>
/// @brief Run the propagation through
/// a uniform energy spectrum and record the
/// energy vs position histograms for each layer
int main() {
    Acts::Logging::Level logLevel = Acts::Logging::VERBOSE;

    // setup the sequencer first w/ config derived from options
    Sequencer::Config seqCfg;
    seqCfg.events = 10;
    seqCfg.numThreads = -1;
    Sequencer sequencer(seqCfg);
    using namespace LUXENavigator;
    Acts::GeometryContext gctx;
    std::string gdmlPath = "lxgeomdump_stave_positron.gdml";
    std::vector<std::string> names{"OPPPSensitive"};
    LUXEGeometry::GeometryOptions gOpt;

    LUXEROOTReader::LUXEROOTSimDataReader::Config readerCfg
        = LUXEROOTReader::defaultSimConfig();
    readerCfg.dataCollection = "SourceLink";
    std::string pathToDir = "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_dataInRootFormat/SignalNextTrial_e1gpc_10.0_1";
    Acts::MagneticFieldContext magCtx;
    // map (x,y,z) -> (x,y,z)
    auto transformPos = [](const Acts::Vector3& pos) {
        return pos;
    };

    // map (Bx,By,Bz) -> (Bx,By,Bz)
    auto transformBField = [](const Acts::Vector3& field, const Acts::Vector3&) {
        return field;
    };

    const std::vector<unsigned int> bins{10u, 10u, 10u};

    auto BField = LUXEMagneticField::buildLUXEBField(transformPos, transformBField, bins);
    std::cout<<BField.getField(Acts::Vector3{3,1,1}).value()<<std::endl;
    auto BFieldPtr = std::make_shared<LUXEMagneticField::BField_t>(BField);

    // Build the LUXE detector
    auto positronArmBpr = LUXEGeometry::makeBlueprint(gdmlPath, names, gctx, gOpt);

    for (const auto & entry : std::filesystem::directory_iterator(pathToDir)) {
        std::string pathToFile = entry.path();
        readerCfg.filePaths.push_back(pathToFile);
    }

    // The events are not sorted in the directory
    // but we need to process them in order
    std::sort(readerCfg.filePaths.begin(), readerCfg.filePaths.end(),
        [] (const std::string& a, const std::string& b) {
            std::size_t idxRootA = a.find_last_of('.');
            std::size_t idxEventA = a.find_last_of('t', idxRootA);
            std::string eventSubstrA = a.substr(idxEventA + 1, idxRootA - idxEventA);

            std::size_t idxRootB = b.find_last_of('.');
            std::size_t idxEventB = b.find_last_of('t', idxRootB);
            std::string eventSubstrB = b.substr(idxEventB + 1, idxRootB - idxEventB);

            return std::stoul(eventSubstrA) < std::stoul(eventSubstrB);
        }
    );

    readerCfg.filePaths = std::vector<std::string>(
        readerCfg.filePaths.begin(), readerCfg.filePaths.begin() + 72);

    // readerCfg.filePaths = {"/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/ActsLUXEPipeline_dataInRootFormat/SignalNextTrial_e1gpc_10.0_1/dataFile_Signal_e1gpc_10.0_EFieldV10p7p1pyN17Vpercm_Processed_Stave25_Event83.root"};

    sequencer.addReader(
        std::make_shared<LUXEROOTReader::LUXEROOTSimDataReader>(readerCfg, logLevel));
    auto testParams = LUXENavigator::makeParameters();

    MeasurementResolution resPixel = {MeasurementType::eLoc01,
                                      {LUXEGeometry::chipSizeX,
                                       LUXEGeometry::chipSizeY}};

    Acts::ViewConfig pConfig = Acts::s_viewSensitive;

    Acts::ObjVisualization3D volumeObj;
    std::vector<std::pair<Acts::GeometryIdentifier,MeasurementResolution>> m;
    auto detector =
            LUXEGeometry::buildLUXEDetector(std::move(positronArmBpr), gctx, gOpt);

    for (auto& vol : detector->rootVolumes()) {
        std::cout<<"Surfaces size: "<<vol->surfaces().size()<<std::endl;
        Acts::GeometryView3D::drawDetectorVolume(
                volumeObj, *(vol), gctx,
                Acts::Transform3::Identity(), pConfig);
        for (auto& surf : vol->surfaces()) {
            std::cout<<"Assigning resolution to surface ID: "<<surf->geometryId()<<std::endl;
            Acts::GeometryView3D::drawSurface(
                    volumeObj, *(surf), gctx,
                    Acts::Transform3::Identity(), pConfig);
            m.push_back(std::make_pair(surf->geometryId(),resPixel));
            std::cout<<"Surface x transform: "<<surf->center(gctx)[0]<<std::endl;
            std::cout<<"Surface y transform: "<<surf->center(gctx)[1]<<std::endl;
            std::cout<<"Surface bounds: "<<surf->bounds()<<std::endl;
        }
    }



    MeasurementResolutionMap resolutions = m;

    auto propagator = LUXENavigator::makePropagator<Acts::EigenStepper<>>(detector, BFieldPtr);

    auto test = LUXENavigator::createMeasurements(propagator, gctx, magCtx, testParams, resolutions);

//    Acts::GeometryView3D::drawSegment(
//            volumeObj, *(surf), gctx,
//            Acts::Transform3::Identity(), pConfig);
    for (auto& sl:test.sourceLinks) {
        std::cout<<sl.parameters<<std::endl;
    }
    volumeObj.write("volumes.obj");

    IdealSeeder::Config seederCfg;
    // seederCfg.roadWidth = 200;
    seederCfg.inputSourceLinks = "SourceLink";
    sequencer.addAlgorithm(
        std::make_shared<IdealSeeder>(seederCfg, logLevel));

    // Run all configured algorithms and return the appropriate status.
    return sequencer.run();
}