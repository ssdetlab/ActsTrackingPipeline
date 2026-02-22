// AlignmentRun.hpp
#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <string>
#include <vector>
#include <toml.hpp>

namespace TrackingPipeline {

struct AlignmentRunConfig {
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

  struct AlignmentStoreInput {
    bool        useFile  = false;
    std::string filePath;
    std::string treeName = "alignment-parameters";
    double      shiftY_mm = 0.0;
    double      shiftZ_mm = 0.0;
  } alignmentStoreInput;

  struct Reader {
    std::vector<std::string> inputFiles;
    std::string treeName           = "fitted-tracks";
    std::string outputMeasurements = "SimMeasurements";
    std::string outputSeedsGuess   = "SeedsGuess";
    std::string outputSeedsEst     = "SeedsEst";
    double      minChi2            = 0.0;
    double      maxChi2            = 18.0;
    bool        mergeIntoOneEvent  = true;
  } reader;

  struct AlignmentKF {
    std::size_t maxSteps        = 1000;

    // origin prior (Acts units, set in .cpp)
    double originStdDevLoc0;
    double originStdDevLoc1;
    double originStdDevPhi;
    double originStdDevTheta;
    double originStdDevQOverP;
    double originStdDevTime;
  } alignmentKF;

  struct AlignmentAlgo {
    std::string inputTrackCandidates      = "SeedsEst";
    std::string outputAlignmentParameters = "AlignmentParameters";
    double      chi2ONdfCutOff            = 1e-16;
    double      deltaChi2ONdfCutOff0      = 50.0;
    double      deltaChi2ONdfCutOff1      = 1e-5;
    std::size_t maxNumIterations          = 200;
    std::string alignmentMode             = "global";
    std::string propagationDirection      = "forward";
  } alignmentAlgo;

  struct Fitting {
    std::size_t maxSteps        = 100000;
    bool        resolvePassive  = false;
    bool        resolveMaterial = true;
    bool        resolveSensitive = true;
    std::string inputCandidates = "SeedsGuess";
    std::string outputTracks    = "Tracks";
  } fitting;

  struct SeedWriter {
    bool        enable             = true;
    std::string inputSeeds         = "SeedsGuess";
    std::string inputTruthClusters = "SimClusters";
    std::string treeName           = "seeds";
    std::string outputFile         = "seeds.root";
  } seedWriter;

  struct TrackWriter {
    bool        enable     = true;
    std::string inputTracks = "Tracks";
    std::string treeName   = "fitted-tracks";
    std::string outputFile = "fitted-tracks.root";
  } trackWriter;

  struct AlignmentWriter {
    bool        enable                    = true;
    std::string inputAlignmentParameters  = "AlignmentParameters";
    std::string treeName                  = "alignment-parameters";
    std::string outputFile                = "alignment-parameters.root";
  } alignmentWriter;

  struct Constraints {
    int minSensitiveIdForConstraints = 40;
    int minSensitiveIdForAlign       = 10;
    int maxSensitiveIdForAlign       = 40;
  } constraints;
};

AlignmentRunConfig parseAlignmentRunConfig(const std::string& path);

Acts::Logging::Level getLogLevel(const std::string& levelStr);

int runAlignment(const std::string& configPath);

} // namespace TrackingPipeline
