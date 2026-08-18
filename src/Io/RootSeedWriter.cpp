#include "TrackingPipeline/Io/RootSeedWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <vector>

#include "TLorentzVector.h"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

RootSeedWriter::RootSeedWriter(const Config& config, Acts::Logging::Level level)
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
  int bufSize = 32000;
  int splitLvl = 0;

  // Measurements
  m_tree->Branch("measurementsGlob", &m_seedMeasurementsGlob, bufSize,
                 splitLvl);
  m_tree->Branch("measurementsLoc", &m_seedMeasurementsLoc, bufSize, splitLvl);
  m_tree->Branch("geoIds", &m_geoIds, bufSize, splitLvl);

  // Seed properties
  m_tree->Branch("eventId", &m_eventId, bufSize, splitLvl);
  m_tree->Branch("size", &m_size, bufSize, splitLvl);
  m_tree->Branch("originMomentumEst", &m_originMomentumEst, bufSize, splitLvl);
  m_tree->Branch("vertexEst", &m_vertexEst, bufSize, splitLvl);
  m_tree->Branch("trackId", &m_trackId, bufSize, splitLvl);

  //------------------------------------------------------------------
  // Initialize the data handles
  m_inputSeeds.initialize(m_cfg.inputSeeds);
  m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
}

ProcessCode RootSeedWriter::finalize() {
  if (m_file != nullptr) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootSeedWriter::write(const AlgorithmContext& ctx) {
  const auto& inputSeeds = m_inputSeeds(ctx);
  const auto& inputTrackParameters = m_inputTrackParameters(ctx);
  const auto& inputSourceLinks = m_inputSourceLinks(ctx);

  ACTS_DEBUG("Received " << inputSeeds.size() << " input seeds");
  ACTS_DEBUG("Received " << inputTrackParameters.size()
                         << " input track parameters");
  ACTS_DEBUG("Received " << inputSourceLinks.size() << " input source links");

  if (inputSeeds.empty()) {
    ACTS_DEBUG("Received empty seed vector. Skipping");
    return ProcessCode::SUCCESS;
  }

  std::lock_guard<std::mutex> lock(m_mutex);

  m_eventId = ctx.eventNumber;

  for (const auto& seed : inputSeeds) {
    const auto& sourceLinkIndices = seed.sourceLinkIndices;
    m_size = sourceLinkIndices.size();

    m_seedMeasurementsGlob.clear();
    m_seedMeasurementsLoc.clear();
    m_geoIds.clear();

    m_seedMeasurementsGlob.reserve(m_size);
    m_seedMeasurementsLoc.reserve(m_size);
    m_geoIds.reserve(m_size);
    for (std::size_t idx : sourceLinkIndices) {
      const auto& ssl = inputSourceLinks.at(idx).get<SimpleSourceLink>();

      const Acts::Vector3& parsGlob = ssl.parametersGlob();
      m_seedMeasurementsGlob.emplace_back(parsGlob.x(), parsGlob.y(),
                                          parsGlob.z());

      const Acts::Vector2& parsLoc = ssl.parametersLoc();
      m_seedMeasurementsLoc.emplace_back(parsLoc.x(), parsLoc.y());
      m_geoIds.push_back(ssl.geometryId().sensitive());
    }

    const auto& originParameters =
        inputTrackParameters.at(seed.originParametersIndex);
    m_originMomentumEst.SetPxPyPzE(
        originParameters.momentum().x(), originParameters.momentum().y(),
        originParameters.momentum().z(), originParameters.absoluteMomentum());
    m_vertexEst.SetXYZ(originParameters.position().x(),
                       originParameters.position().y(),
                       originParameters.position().z());

    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
