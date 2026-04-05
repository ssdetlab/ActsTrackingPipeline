#include "TrackingPipeline/Io/RootMeasurementWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include <Acts/Utilities/Logger.hpp>

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

RootMeasurementWriter::RootMeasurementWriter(const Config& config,
                                             Acts::Logging::Level level)
    : m_cfg(config), m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_file = new TFile(m_cfg.filePath.c_str(), "RECREATE");
  m_tree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());

  //------------------------------------------------------------------
  // Tree branches
  int buf_size = 32000;
  int split_lvl = 0;

  // Parameters at measurements
  m_tree->Branch("geoCenterGlobal", &m_geoCenterGlobal, buf_size, split_lvl);
  m_tree->Branch("geoCenterLocal", &m_geoCenterLocal, buf_size, split_lvl);
  m_tree->Branch("cov", &m_cov, buf_size, split_lvl);
  m_tree->Branch("geoId", &m_geoId, buf_size, split_lvl);
  m_tree->Branch("eventId", &m_eventId, buf_size, split_lvl);

  //------------------------------------------------------------------
  // Initialize the data handles
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_inputSourceLinkIndices.initialize(m_cfg.inputSourceLinkIndices);
}

ProcessCode RootMeasurementWriter::finalize() {
  if (m_file != nullptr) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootMeasurementWriter::write(const AlgorithmContext& ctx) {
  const auto& inputSourceLinks = m_inputSourceLinks(ctx);
  const auto& inputSourceLinkIndices = m_inputSourceLinkIndices(ctx);

  ACTS_DEBUG("Received " << inputSourceLinks.size() << " source links");
  ACTS_DEBUG("Received " << inputSourceLinkIndices.size()
                         << " source link indices");

  std::lock_guard<std::mutex> lock(m_mutex);

  for (std::size_t i = 0; i < inputSourceLinkIndices.size(); i++) {
    std::size_t idx = inputSourceLinkIndices.at(i);
    const auto& ssl = inputSourceLinks.at(idx).get<SimpleSourceLink>();

    Acts::Vector2 geoCenterLocal = ssl.parametersLoc();
    Acts::Vector3 geoCenterGlobal = ssl.parametersGlob();

    m_geoCenterGlobal.SetX(geoCenterGlobal.x());
    m_geoCenterGlobal.SetY(geoCenterGlobal.y());
    m_geoCenterGlobal.SetZ(geoCenterGlobal.z());

    m_geoCenterLocal.SetX(geoCenterLocal.x());
    m_geoCenterLocal.SetY(geoCenterLocal.y());

    m_geoId = ssl.geometryId().sensitive();
    m_eventId = ssl.eventId();

    m_cov(0, 0) = ssl.covariance()(0, 0);
    m_cov(0, 1) = ssl.covariance()(0, 1);
    m_cov(1, 0) = ssl.covariance()(1, 0);
    m_cov(1, 1) = ssl.covariance()(1, 1);

    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
