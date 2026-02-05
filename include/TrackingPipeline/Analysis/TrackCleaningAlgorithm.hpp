#pragma once

#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp"

#include <string>

class TrackCleaningAlgorithm final : public IAlgorithm {
 public:
  using CleaningContainer = CleaningTracks; // std::vector<TrackDescriptor>
  using ActsContainer     = Tracks;         // wrapper around ActsTracks

  struct Config {
    // Optional: pre-flattened cleaning tracks (from RootTrackReader)
    std::string inputCleaningTracks;
    // Optional: Acts tracks (from KFTrackFittingAlgorithm)
    std::string inputActsTracks;
    // Required: cleaned CleaningTracks collection
    std::string outputTracks;
  };

  TrackCleaningAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const AlgorithmContext& ctx) const override;

 private:
  Config m_cfg;

  ReadDataHandle<CleaningContainer> m_inputCleaning{this, "InputCleaningTracks"};
  ReadDataHandle<ActsContainer>     m_inputActs{this, "InputActsTracks"};
  WriteDataHandle<CleaningContainer> m_outputTracks{this, "OutputTracks"};
};

