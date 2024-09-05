#include "ActsLUXEPipeline/Sequencer.hpp"
#include "ActsLUXEPipeline/E320Geometry.hpp"
#include "ActsLUXEPipeline/BinnedMagneticField.hpp"
#include "ActsLUXEPipeline/DipoleMagField.hpp"
#include "ActsLUXEPipeline/QuadrupoleMagField.hpp"
#include "ActsLUXEPipeline/CompositeMagField.hpp"
#include "ActsLUXEPipeline/MeasurementsCreator.hpp"
#include "ActsLUXEPipeline/AlgorithmContext.hpp"
#include "ActsLUXEPipeline/ROOTLookupDataWriter.hpp"
#include "ActsLUXEPipeline/CsvLookupTableWriter.hpp"
#include "ActsLUXEPipeline/Generators.hpp"

#include "Acts/Utilities/Logger.hpp"

#include <filesystem>

using namespace Acts::UnitLiterals;

using Propagator = Acts::Propagator<
    Acts::EigenStepper<>, 
    Acts::Experimental::DetectorNavigator>;
using TrackParameters = Acts::CurvilinearTrackParameters;

/// @brief Run the propagation through 
/// a uniform energy spectrum and record the
/// energy vs position histograms for each layer
int main() {
    // Set the log level
    Acts::Logging::Level logLevel = Acts::Logging::INFO;

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
        "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/E320Pipeline_gdmls/ettgeom_magnet_pdc_tracker.gdml";
    std::vector<std::string> names{"OPPPSensitive", "DetChamberWindow"};

    std::string materialPath = "/home/romanurmanov/lab/LUXE/acts_LUXE_tracking/E320Pipeline_material/uniform/material.json";

    // Build the detector
    auto trackerBP = 
        E320Geometry::makeBlueprintE320(gdmlPath, names, gOpt);
    auto detector =
        E320Geometry::buildE320Detector(std::move(trackerBP), gctx, gOpt, {});

    // --------------------------------------------------------------
    // The magnetic field setup

    // Extent in already rotated frame
    Acts::Extent quad1Extent;
    quad1Extent.set(
        Acts::BinningValue::binX, 
        gOpt.quad1Translation[0] - gOpt.quad1Bounds[0],
        gOpt.quad1Translation[0] + gOpt.quad1Bounds[0]);
    quad1Extent.set(
        Acts::BinningValue::binZ,
        gOpt.quad1Translation[1] - gOpt.quad1Bounds[1],
        gOpt.quad1Translation[1] + gOpt.quad1Bounds[1]);
    quad1Extent.set(
        Acts::BinningValue::binY,
        gOpt.quad1Translation[2] - gOpt.quad1Bounds[2],
        gOpt.quad1Translation[2] + gOpt.quad1Bounds[2]);

    Acts::Extent quad2Extent;
    quad2Extent.set(
        Acts::BinningValue::binX, 
        gOpt.quad2Translation[0] - gOpt.quad2Bounds[0],
        gOpt.quad2Translation[0] + gOpt.quad2Bounds[0]);
    quad2Extent.set(
        Acts::BinningValue::binZ,
        gOpt.quad2Translation[1] - gOpt.quad2Bounds[1],
        gOpt.quad2Translation[1] + gOpt.quad2Bounds[1]);
    quad2Extent.set(
        Acts::BinningValue::binY,
        gOpt.quad2Translation[2] - gOpt.quad2Bounds[2],
        gOpt.quad2Translation[2] + gOpt.quad2Bounds[2]);

    Acts::Extent quad3Extent;
    quad3Extent.set(
        Acts::BinningValue::binX, 
        gOpt.quad3Translation[0] - gOpt.quad3Bounds[0],
        gOpt.quad3Translation[0] + gOpt.quad3Bounds[0]);
    quad3Extent.set(
        Acts::BinningValue::binZ,
        gOpt.quad3Translation[1] - gOpt.quad3Bounds[1],
        gOpt.quad3Translation[1] + gOpt.quad3Bounds[1]);
    quad3Extent.set(
        Acts::BinningValue::binY,
        gOpt.quad3Translation[2] - gOpt.quad3Bounds[2],
        gOpt.quad3Translation[2] + gOpt.quad3Bounds[2]);

    Acts::Extent dipoleExtent;
    dipoleExtent.set(
        Acts::BinningValue::binX, 
        gOpt.dipoleTranslation.x() - gOpt.dipoleBounds[0],
        gOpt.dipoleTranslation.x() + gOpt.dipoleBounds[0]);
    dipoleExtent.set(
        Acts::BinningValue::binZ,
        gOpt.dipoleTranslation.y() - gOpt.dipoleBounds[1],
        gOpt.dipoleTranslation.y() + gOpt.dipoleBounds[1]);
    dipoleExtent.set(
        Acts::BinningValue::binY,
        gOpt.dipoleTranslation.z() - gOpt.dipoleBounds[2],
        gOpt.dipoleTranslation.z() + gOpt.dipoleBounds[2]);

    QuadrupoleMagField quad1Field(
        gOpt.quadrupolesParams[0], 
        gOpt.actsToWorldRotation.inverse() * gOpt.quad1Translation, 
        gOpt.actsToWorldRotation);
    QuadrupoleMagField quad2Field(
        gOpt.quadrupolesParams[1], 
        gOpt.actsToWorldRotation.inverse() * gOpt.quad2Translation, 
        gOpt.actsToWorldRotation);
    QuadrupoleMagField quad3Field(
        gOpt.quadrupolesParams[2], 
        gOpt.actsToWorldRotation.inverse() * gOpt.quad3Translation, 
        gOpt.actsToWorldRotation);
    
    Acts::ActsScalar dipoleB = 0.31_T;
    DipoleMagField dipoleField(
        gOpt.dipoleParams, 
        dipoleB,
        gOpt.actsToWorldRotation, 
        gOpt.actsToWorldRotation.inverse() * gOpt.dipoleTranslation);
    
    CompositeMagField::FieldComponents fieldComponents = {
        {quad1Extent, &quad1Field},
        {quad2Extent, &quad2Field},
        {quad3Extent, &quad3Field},
        {dipoleExtent, &dipoleField}
    };

    auto field = std::make_shared<CompositeMagField>(fieldComponents);

    // --------------------------------------------------------------
    // Event creation 

    // Setup the sequencer
    Sequencer::Config seqCfg;
    seqCfg.events = 200000;
    seqCfg.numThreads = 1;
    seqCfg.trackFpes = false;
    Sequencer sequencer(seqCfg);

    // Setup the measurements creator
    Acts::Experimental::DetectorNavigator::Config navCfg;
    navCfg.detector = detector.get();
    navCfg.resolvePassive = false;
    navCfg.resolveMaterial = true;
    navCfg.resolveSensitive = true;

    Acts::Experimental::DetectorNavigator navigator(
        navCfg, 
        Acts::getDefaultLogger(
            "DetectorNavigator", 
            logLevel));
    Acts::EigenStepper<> stepper(field);

    auto propagator = 
        Propagator(std::move(stepper), std::move(navigator));

    RangedUniformMomentumGenerator momGen;
    momGen.Pranges = {
        {0.5_GeV, 1.0_GeV},
        {1.0_GeV, 1.5_GeV},
        {1.5_GeV, 2.0_GeV},
        {2.0_GeV, 2.5_GeV},
        {2.5_GeV, 3.0_GeV},
        {3.0_GeV, 3.5_GeV},
        {3.5_GeV, 4.0_GeV},
        {4.0_GeV, 4.5_GeV}};

    // UniformVertexGenerator vertexGen;
    // vertexGen.mins = {-100_um, -100_um, -100_um};
    // vertexGen.maxs = {100_um, 100_um, 100_um};

    // GaussianMomentumGenerator momGen; 
    // momGen.pMagRange = {0.5_GeV, 4.5_GeV};
    // momGen.thetaRange = {-M_PI / 4, M_PI / 4};
    // momGen.phiRange = {0, 2 * M_PI};

    MeasurementsCreator::Config mcCfg;
    mcCfg.outputCollection = "Measurements";
    mcCfg.vertexGenerator = std::make_shared<StationaryVertexGenerator>();
    mcCfg.momentumGenerator = std::make_shared<RangedUniformMomentumGenerator>(momGen);
    mcCfg.randomNumberSvc = std::make_shared<RandomNumbers>(RandomNumbers::Config());
    mcCfg.nTracks = 1;

    sequencer.addAlgorithm(
        std::make_shared<MeasurementsCreator>(
            propagator, mcCfg, logLevel));

    // --------------------------------------------------------------
    // Lookup data generation 

    // Add the lookup data writers
    ROOTLookupDataWriter::Config lookupWriterCfg;
        
    lookupWriterCfg.inputCollection = mcCfg.outputCollection;

    SimpleSourceLink::SurfaceAccessor surfaceAccessor{*detector};
    lookupWriterCfg.surfaceAccessor.connect<
        &SimpleSourceLink::SurfaceAccessor::operator()>(
        &surfaceAccessor);

    // Extent in already rotated frame
    Acts::Extent firstLayerExtent;
    firstLayerExtent.set(
        Acts::BinningValue::binX, 
        gOpt.chipX - gOpt.chipSizeX/2 - 1_mm,
        gOpt.chipX + gOpt.chipSizeX/2 + 1_mm);
    firstLayerExtent.set(
        Acts::BinningValue::binZ,
        -gOpt.chipY.at(8) - gOpt.chipSizeY/2 - 1_mm,
        -gOpt.chipY.at(0) + gOpt.chipSizeY/2 + 1_mm);
    firstLayerExtent.set(
        Acts::BinningValue::binY,
        gOpt.layerZPositions.at(0) - gOpt.layerBounds.at(2) - 1_mm,
        gOpt.layerZPositions.at(0) + gOpt.layerBounds.at(2) + 1_mm);

    lookupWriterCfg.firstLayerExtent = firstLayerExtent;

    sequencer.addWriter(
        std::make_shared<ROOTLookupDataWriter>(lookupWriterCfg, logLevel));

    auto lookupMakerCfg = CsvLookupTableWriter::Config();

    lookupMakerCfg.filePath = "lookupTable.csv";

    // Grid parameters
    lookupMakerCfg.YFirst = {
        1000, 
        -gOpt.chipY.at(8) - gOpt.chipSizeY/2 - 1_mm, 
        -gOpt.chipY.at(0) + gOpt.chipSizeY/2 + 1_mm};
    lookupMakerCfg.XFirst = {
        100, 
        gOpt.chipX - gOpt.chipSizeX/2 - 1_mm,
        gOpt.chipX + gOpt.chipSizeX/2 + 1_mm};
    lookupMakerCfg.surfaceAccessor.connect<
        &SimpleSourceLink::SurfaceAccessor::operator()>(
        &surfaceAccessor);
    lookupMakerCfg.firstLayerExtent = firstLayerExtent;

    auto lookupMaker = std::make_shared<CsvLookupTableWriter>(lookupMakerCfg, logLevel);

    sequencer.addWriter(lookupMaker);

    return sequencer.run();
}
