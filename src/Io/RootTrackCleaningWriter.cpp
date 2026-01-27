#include "TrackingPipeline/Io/RootTrackCleaningWriter.hpp"

#include "TrackingPipeline/Infrastructure/WriterRegistry.hpp"

#include <TFile.h>
#include <TTree.h>

#include <toml.hpp>

RootCleaningTrackWriter::RootCleaningTrackWriter(
    const Config& cfg, Acts::Logging::Level level)
  : IWriter()
  , m_cfg(cfg)
  , m_logger(Acts::getDefaultLogger(name(), level)) {

  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("RootCleaningTrackWriter: inputTracks empty");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("RootCleaningTrackWriter: treeName empty");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("RootCleaningTrackWriter: filePath empty");
  }

  m_inputTracks.initialize(m_cfg.inputTracks);
}

void RootCleaningTrackWriter::initializeTree() {
  if (m_fileInitialized) {
    return;
  }

  m_file = TFile::Open(m_cfg.filePath.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie()) {
    throw std::runtime_error("RootCleaningTrackWriter: cannot create file " +
                             m_cfg.filePath);
  }

  m_tree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());

  // Simple branches: basic track info and hit coordinates
  m_tree->Branch("eventId", &m_eventId);
  m_tree->Branch("trackId", &m_trackId);
  m_tree->Branch("pdgId",   &m_pdgId);
  m_tree->Branch("charge",  &m_charge);
  m_tree->Branch("chi2Smoothed", &m_chi2);
  m_tree->Branch("ndf",     &m_ndf);

  m_tree->Branch("hitX", &m_hitX);
  m_tree->Branch("hitY", &m_hitY);
  m_tree->Branch("hitZ", &m_hitZ);

  m_fileInitialized = true;
}

ProcessCode RootCleaningTrackWriter::write(const AlgorithmContext& ctx) {
  const auto& tracks = m_inputTracks(ctx);

  if (!m_fileInitialized) {
    initializeTree();
  }

  // Clear buffers
  m_hitX.clear();
  m_hitY.clear();
  m_hitZ.clear();

  for (const auto& trk : tracks) {
    m_eventId = trk.eventId;
    m_trackId = trk.trackId;
    m_pdgId   = trk.pdgId;
    m_charge  = trk.charge;
    m_chi2    = trk.chi2Smoothed;
    m_ndf     = trk.ndf;

    m_hitX.clear();
    m_hitY.clear();
    m_hitZ.clear();
    m_hitX.reserve(trk.trackHitsGlobal.size());
    m_hitY.reserve(trk.trackHitsGlobal.size());
    m_hitZ.reserve(trk.trackHitsGlobal.size());

    for (const auto& pos : trk.trackHitsGlobal) {
      m_hitX.push_back(pos.x());
      m_hitY.push_back(pos.y());
      m_hitZ.push_back(pos.z());
    }

    m_tree->Fill();
  }

  return ProcessCode::SUCCESS;
}

RootCleaningTrackWriter::~RootCleaningTrackWriter() {
  if (m_fileInitialized && m_file) {
    m_file->cd();
    if (m_tree) {
      m_tree->Write("", TObject::kOverwrite);
    }
    m_file->Close();
    m_file = nullptr;
    m_tree = nullptr;
    m_fileInitialized = false;
  }
}



namespace {
struct RootCleaningTrackWriterRegistrar {
  RootCleaningTrackWriterRegistrar() {
    using namespace TrackingPipeline;

    WriterRegistry::instance().registerBuilder(
        "RootCleaningTrackWriter",
        [](const toml::value& section,
           Acts::Logging::Level logLevel,
           const std::string& runRoot) -> WriterPtr {

          RootCleaningTrackWriter::Config cfg;
          cfg.inputTracks =
              toml::find<std::string>(section, "inputTracks");
          cfg.treeName =
              toml::find<std::string>(section, "treeName");
          auto relPath =
              toml::find<std::string>(section, "filePath");
          cfg.filePath = runRoot + "/" + relPath;

          return std::make_shared<RootCleaningTrackWriter>(cfg, logLevel);
        });
  }
} _RootCleaningTrackWriterRegistrar;
}  // namespace
