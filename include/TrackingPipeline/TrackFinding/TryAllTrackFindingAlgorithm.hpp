#pragma once

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/detail/SeedTree.hpp"

class TryAllTrackFindingAlgorithm : public IAlgorithm {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// The input collection
    std::string inputCollection = "Seed";
    /// The output collection
    std::string outputCollection = "Seed";
    /// Minimum number of source links
    int minSourceLinks = 3;
    /// Maximum number of source links
    int maxSourceLinks = 10;
  };

  /// @brief Constructor
  TryAllTrackFindingAlgorithm(Config config, Acts::Logging::Level level)
      : IAlgorithm("TryAllTrackFindingAlgorithm", level),
        m_cfg(std::move(config)) {
    m_inputSeeds.initialize(m_cfg.inputCollection);
    m_outputSeeds.initialize(m_cfg.outputCollection);
  }
  ~TryAllTrackFindingAlgorithm() = default;

  /// @brief The execute method
  ProcessCode execute(const AlgorithmContext& ctx) const override {
    // Get the input seeds
    // from the context
    auto input = m_inputSeeds(ctx);

    std::cout << "TRACK FINDING ALL: INPUT SEEDS SIZE = " << input.size()
              << std::endl;

    if (input.empty()) {
      m_outputSeeds(ctx, Seeds());
      return ProcessCode::SUCCESS;
    }

    Seeds trackCandidates;

    for (const auto& seed : input) {
      if (seed.sourceLinks.size() < m_cfg.minSourceLinks ||
          seed.sourceLinks.size() > m_cfg.maxSourceLinks) {
        continue;
      }
      std::vector<std::vector<Acts::SourceLink>> tracks;
      SeedTree seedTree(seed);
      constructTracks(seedTree.root, {}, tracks);

      for (const auto& track : tracks) {
        trackCandidates.push_back(Seed{track, seed.ipParameters, seed.trackId});
      }
    }

    std::cout << "TRACK FINDING ALL: OUT SEEDS SIZE = "
              << trackCandidates.size() << std::endl;
    m_outputSeeds(ctx, std::move(trackCandidates));

    return ProcessCode::SUCCESS;
  }

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<Seeds> m_inputSeeds{this, "InputSeeds"};

  WriteDataHandle<Seeds> m_outputSeeds{this, "OutputSeeds"};

  void constructTracks(
      std::shared_ptr<SeedTree::Node> root, std::vector<Acts::SourceLink> track,
      std::vector<std::vector<Acts::SourceLink>>& tracks) const {
    if (root->children.size() == 0) {
      track.push_back(root->m_sourceLink);
      tracks.push_back(track);
    }

    track.push_back(root->m_sourceLink);
    for (auto& child : root->children) {
      constructTracks(child, track, tracks);
    }
  }
};
