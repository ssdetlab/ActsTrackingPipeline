#include "TrackingPipeline/Io/RootMeasurementWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"

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
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
}

ProcessCode RootMeasurementWriter::finalize() {
  if (m_file) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootMeasurementWriter::write(const AlgorithmContext& ctx) {
  auto inputMeasurements = m_inputMeasurements(ctx);

  std::lock_guard<std::mutex> lock(m_mutex);

  for (const auto& meas : inputMeasurements) {
    const auto& ssl = meas.get<SimpleSourceLink>();
    const auto* surf = m_cfg.surfaceAccessor(meas);

    Acts::Vector2 geoCenterLocal = ssl.parametersLoc();
    Acts::Vector3 geoCenterGlobal = surf->localToGlobal(
        ctx.geoContext, ssl.parametersLoc(), Acts::Vector3(0, 1, 0));

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
