// include/TrackingPipeline/Run/FastSimRun.hpp
#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <string>
#include <toml.hpp>

namespace TrackingPipeline {

struct FastSimRunConfig {
  toml::value root;

  struct Main {
    std::string logLevel = "INFO";
    std::string outputDir;
  } main;

  struct Sequencer {
    std::size_t numThreads = 1;
    std::size_t skip       = 0;
    bool        trackFpes  = false;
  } sequencer;

  struct Alignment {
    bool        useFile  = false;
    std::string filePath;
    std::string treeName = "alignment-parameters";
  } alignment;

  struct DummyReader {
    std::size_t nEvents           = 100000;
    std::string outputSourceLinks = "SimMeasurements";
    std::string outputSimClusters = "SimClusters";
  } dummyReader;

  struct Signal {
    bool resolvePassive        = false;
    bool resolveMaterial       = true;
    bool resolveSensitive      = true;
    double resY_um     = 5.0;   // digitizer
    double resZ_um     = 5.0;   // digitizer
    double pMin_GeV    = 2.0;
    double pMax_GeV    = 3.0;
    double phiMin      = -0.001;
    double phiMax      =  0.001;
    double thetaMin    = 1.57079632679 - 0.003;
    double thetaMax    = 1.57079632679 + 0.003;
    double vtxSigma_um = 30.0;
  } signal;

  struct SignalMeasCreator {
    std::string inputSourceLinks   = "SimMeasurements";
    std::string inputSimClusters   = "SimClusters";
    std::string outputSourceLinks  = "Measurements";
    std::string outputSimClusters  = "Clusters";
    std::size_t maxSteps           = 1000;
    int         nMeasurements      = 1;
    int         minSensitiveId     = 40;         // constraints gid.sensitive() >= this
  } signalMeasCreator;

  struct Background {
    bool        enable   = false;
    std::size_t nHits    = 700;
    double      resY_um  = 5.0;
    double      resZ_um  = 5.0;

    std::string inputSourceLinks  = "Measurements1";
    std::string inputSimClusters  = "Clusters1";
    std::string outputSourceLinks = "Measurements";
    std::string outputSimClusters = "Clusters";
    int         nMeasurements     = 1;
  } background;

  struct ClusterWriter {
    bool        enable     = true;
    std::string inputClusters = "Clusters";
    std::string treeName   = "clusters";
    std::string filePrefix = "clusters-";
  } clusterWriter;
};

FastSimRunConfig parseFastSimRunConfig(const std::string& path);

Acts::Logging::Level getLogLevel(const std::string& levelStr);

int runFastSim(const std::string& configPath, int jobId);

} // namespace TrackingPipeline
