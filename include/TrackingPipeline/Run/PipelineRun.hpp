#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>
#include <string>
#include <vector>
#include <toml.hpp>

namespace TrackingPipeline {

using SurfaceMap = std::map<Acts::GeometryIdentifier, const Acts::Surface*>;

struct PipelineConfig {
  toml::value root;

  struct Main {
    std::string outputDirLoc;
    std::string outputDirName;
    std::string logLevel;
  } main;

  struct Detector {
    bool fillSurfaceMap;
    // Alignment store configuration
    std::string storeMode;        // "fixed" or "gaussian"
    std::string longBinningName;  // e.g. "binY"
    std::string shortBinningName; // e.g. "binZ"
    double      longTransStdUm;   // [um]
    double      shortTransStdUm;  // [um]
  } detector;

  struct AlignmentProvider {
    bool        enable = false;
    std::string treeName;
    std::string parametersFilepath;
  } alignmentProvider;

  struct DataReader {
    std::string type;
  } dataReader;

  struct DataWriter {
    std::vector<std::string> types;
  } dataWriter;

  struct Seeding {
    bool enable;
    std::vector<std::string> types;
  } seeding;

  struct TrackFitting {
    bool enable;
    std::vector<std::string> types;
  } trackFitting;

  struct TrackCleaning {
  bool enable;
  std::vector<std::string> types;
  } trackCleaning;

  struct Preprocessing {
    bool enable = false;
    std::vector<std::string> inputDirs;
    std::string inputTreeName;
    std::string inputBranchName = "event";
    std::string outputFile;
    std::string outputTreeName = "MyTree";
    std::size_t skipEntries = 0;
  } preprocessing;

};

PipelineConfig parsePipelineConfig(const std::string& path);

Acts::Logging::Level getLogLevel(const std::string& levelStr);

// main body: returns 0/!=0
int runPipeline(const std::string& configPath);

} // namespace TrackingPipeline
