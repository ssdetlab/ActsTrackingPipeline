#pragma once

#include "Acts/Utilities/Logger.hpp"

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

#include <memory>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

class RootCleaningTrackWriter : public IWriter {
 public:
  struct Config {
    std::string inputTracks;  // e.g. "CleaningTracksFiltered"
    std::string treeName;     // e.g. "fitted-tracks"
    std::string filePath;     // e.g. "tracks-filtered.root"
  };

  RootCleaningTrackWriter(const Config& cfg, Acts::Logging::Level level);

  ~RootCleaningTrackWriter() override;

  std::string name() const override { return "RootCleaningTrackWriter"; }

  ProcessCode write(const AlgorithmContext& ctx) override;

  const Config& config() const { return m_cfg; }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;

  // Reads CleaningTracks from the whiteboard
  ReadDataHandle<CleaningTracks> m_inputTracks{this, "Tracks"};

  std::unique_ptr<const Acts::Logger> m_logger;

  // ROOT objects
  bool m_fileInitialized = false;
  TFile* m_file = nullptr;
  TTree* m_tree = nullptr;

  // Branch buffers
  std::size_t m_eventId = 0;
  std::size_t m_trackId = 0;
  int         m_pdgId   = 0;
  int         m_charge  = 0;
  double      m_chi2    = 0.0;
  std::size_t m_ndf     = 0;

  std::vector<double> m_hitX;
  std::vector<double> m_hitY;
  std::vector<double> m_hitZ;

  void initializeTree();
};
