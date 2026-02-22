// include/TrackingPipeline/Run/FullSimTrackingRun.hpp
#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <string>
#include <toml.hpp>

namespace TrackingPipeline {

struct FullSimTrackingRunConfig {
  toml::value root;

  struct Main {
    std::string logLevel = "INFO";
    std::string outputDir;
  } main;

  struct Sequencer {
    std::size_t numThreads = 1;
    std::size_t skip       = 0;
    bool trackFpes         = false;
  } sequencer;

  struct Reader {
    std::string clustersDir;
    std::string treeName           = "clusters";
    int minGeoId                   = 10;
    int maxGeoId                   = 18;
    bool surfaceLocalToGlobal      = true;
    std::string outputSourceLinks  = "Measurements";
    std::string outputSimClusters  = "SimClusters";
  } reader;

  struct Alignment {
    // If true, read alignment from ROOT; if false, use makeAlignmentStore
    bool        useFile  = false;
    std::string filePath;
    std::string treeName = "alignment-parameters";
  } alignment;

  struct Seeding {
    std::string inputSourceLinks = "Measurements";
    std::string outputSeeds      = "Seeds";
    int         minLayers        = 5;
    int         maxLayers        = 5;
    double      beamlineTilt     = 0.0;

    // Origin uncertainty (Acts units, set in .cpp)
    double originStdDevLoc0;
    double originStdDevLoc1;
    double originStdDevPhi;
    double originStdDevTheta;
    double originStdDevQOverP;
    double originStdDevTime;

    // multiple scattering
    double X0  = 21.82;
    double rho = 2.329;
    double x   = 25e-4;
    double P   = 2500.;
    double z   = 1.;

    // HT grid / options
    int    nCellsThetaShort = 500;
    int    nCellsRhoShort   = 4000;
    int    nCellsThetaLong  = 500;
    int    nCellsRhoLong    = 4000;
    int    nGLSIterations   = 2;
    int    minXCount        = 3;
    int    minSeedSize      = 5;
    int    maxSeedSize      = 100;
    double maxChi2          = 1e16;
  } seeding;

  struct Fitting {
    std::size_t maxSteps        = 100000;
    bool        resolvePassive  = false;
    bool        resolveMaterial = true;
    bool        resolveSensitive = true;
    std::string inputCandidates = "Seeds";
    std::string outputTracks    = "Tracks";
  } fitting;

  struct MeasurementWriter {
    bool        enable     = false;
    std::string input      = "Measurements";
    std::string treeName   = "measurements";
    std::string outputFile = "measurements.root";
  } measurementWriter;

  struct SeedWriter {
    bool        enable             = false;
    std::string inputSeeds         = "Seeds";
    std::string inputTruthClusters = "SimClusters";
    std::string treeName           = "seeds";
    std::string outputFile         = "seeds-0.root";
  } seedWriter;

  struct TrackWriter {
    bool        enable           = true;
    std::string inputTracks      = "Tracks";
    std::string inputSimClusters = "SimClusters";
    std::string treeName         = "fitted-tracks";
    std::string outputFile       = "fitted-tracks-0.root";
  } trackWriter;
};

FullSimTrackingRunConfig parseFullSimTrackingRunConfig(const std::string& path);

Acts::Logging::Level getLogLevel(const std::string& levelStr);

int runFullSimTracking(const std::string& configPath);

} // namespace TrackingPipeline
