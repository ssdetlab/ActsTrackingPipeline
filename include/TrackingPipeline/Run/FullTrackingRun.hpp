// FullTrackingRun.hpp
#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <string>
#include <toml.hpp>

namespace TrackingPipeline {

struct FullTrackingConfig {
  toml::value root;

  struct Main {
    std::string logLevel = "INFO";
    std::string outputDir;
  } main;

  struct Sequencer {
    std::size_t skip       = 0;
    std::size_t numThreads = 1;
    bool trackFpes         = false;
  } sequencer;

  struct Reader {
    std::string inputDir;
    std::string treeName    = "MyTree";
    std::string eventKey    = "event";
    std::string outputSourceLinks = "Measurements";
  } reader;

  struct Seeding {
    std::string inputSourceLinks = "Measurements";
    std::string outputSeeds = "Seeds";
    int minLayers           = 5;
    int maxLayers           = 5;
    double beamlineTilt     = 0.0;

    // Origin uncertainty for seeding
    double originStdDevLoc0; //   = 10.0 * 1_mm;
    double originStdDevLoc1; //   = 10.0 * 1_mm;
    double originStdDevPhi; //    = 10.0 * 1_degree;
    double originStdDevTheta; //  = 10.0 * 1_degree;
    double originStdDevQOverP; // = 0.01 / 1_GeV;
    double originStdDevTime; //   = 1.0 * 1_fs;
    // multiple scattering
    double X0  = 21.82;
    double rho = 2.329;
    double x   = 25e-4;
    double P   = 2500.;
    double z   = 1.;

    // HT grid / options
    int nCellsThetaShort = 500;
    int nCellsRhoShort   = 4000;
    int nCellsThetaLong  = 500;
    int nCellsRhoLong    = 4000;
    int nGLSIterations   = 2;
    int minXCount        = 4;
    int minSeedSize      = 5;
    int maxSeedSize      = 100;
    double maxChi2       = 1e2;
  } seeding;

  struct Fitting {
    std::size_t maxSteps       = 100000;
    bool resolvePassive        = false;
    bool resolveMaterial       = true;
    bool resolveSensitive      = true;
    std::string inputCandidates = "Seeds";
    std::string outputTracks    = "Tracks";
  } fitting;

  struct Alignment {
    bool enable = false;
    std::string filePath;
    std::string treeName = "alignment-parameters";
  } alignment;

  struct MeasurementWriter {
    bool        enable      = true;
    std::string input       = "Measurements";
    std::string treeName    = "measurements";
    std::string outputFile  = "measurements.root";
  } measurementWriter;

  struct SeedWriter {
    bool        enable      = true;
    std::string input       = "Seeds";
    std::string treeName    = "seeds";
    std::string outputFile  = "seeds.root";
  } seedWriter;

  struct TrackWriter {
    bool        enable      = true;
    std::string input       = "Tracks";
    std::string treeName    = "fitted-tracks";
    std::string outputFile  = "fitted-tracks.root";
  } trackWriter;
};

FullTrackingConfig parseFullTrackingConfig(const std::string& path);

Acts::Logging::Level getLogLevel(const std::string& levelStr);

int runFullTracking(const std::string& configPath);

} // namespace TrackingPipeline
