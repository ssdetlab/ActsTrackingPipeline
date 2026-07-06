#include "TrackingPipeline/Io/E320RootDataReader.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"

#include <cstddef>
#include <stdexcept>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

E320::E320RootDataReader::E320RootDataReader(const Config& config,
                                             Acts::Logging::Level level)
    : IReader(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.filePaths.empty()) {
    throw std::invalid_argument("Missing input filenames");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  if (m_cfg.filePaths.size() == 1) {
    m_file = new TFile(m_cfg.filePaths.at(0).c_str());
    m_tree = m_file->Get<TTree>(m_cfg.treeName.c_str());
  } else {
    m_chainOwner = new TChain(m_cfg.treeName.c_str());
    // Add the files to the chain
    for (const auto& path : m_cfg.filePaths) {
      m_chainOwner->Add(path.c_str());
    }
    m_tree = dynamic_cast<TTree*>(m_chainOwner);
  }

  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
  m_outputDetSourceLinksIndices.initialize(m_cfg.outputDetSourceLinkIndices);
  m_outputBpmSourceLinksIndices.initialize(m_cfg.outputBpmSourceLinkIndices);
  m_outputEventMetaData.initialize(m_cfg.outputEventMetaData);

  // Set eventId branch
  m_tree->SetBranchAddress("eventId", &m_eventId);
  if (m_tree->GetBranch("eventId") == nullptr) {
    throw std::invalid_argument("Missing eventId branch");
  }
  m_tree->SetBranchStatus("eventId", true);
  auto nEntries = static_cast<std::size_t>(m_tree->GetEntries());

  // Add the first entry
  m_tree->GetEntry(0);
  m_eventMap.emplace_back(m_eventId, 0, 0);

  // Go through all entries and store the position of the events
  for (std::size_t i = 1; i < nEntries; ++i) {
    m_tree->GetEntry(i);
    if (m_eventId != std::get<0>(m_eventMap.back())) {
      std::get<2>(m_eventMap.back()) = i;
      m_eventMap.emplace_back(m_eventId, i, i);
    }
    if (i == nEntries - 1) {
      std::get<2>(m_eventMap.back()) = nEntries;
    }
  }

  // Sort by event id
  std::ranges::sort(m_eventMap, [](const auto& a, const auto& b) {
    return std::get<0>(a) < std::get<0>(b);
  });

  // Set data branch
  m_tree->SetBranchAddress(m_cfg.eventKey.c_str(), &m_detEvent);

  // Initialize geometry ids
  const auto& goInst = *E320::GeometryOptions::instance();
  m_geoIdMap = {{0, goInst.tcParameters.at(4).geoId},
                {2, goInst.tcParameters.at(3).geoId},
                {4, goInst.tcParameters.at(2).geoId},
                {6, goInst.tcParameters.at(1).geoId},
                {8, goInst.tcParameters.at(0).geoId},
                {3156, goInst.bpm0Parameters.geoId},
                {3218, goInst.bpm1Parameters.geoId},
                {3265, goInst.bpm2Parameters.geoId},
                {3315, goInst.bpm3Parameters.geoId}};

  ACTS_DEBUG("Event range: " << availableEvents().first << " - "
                             << availableEvents().second);
}

std::pair<std::size_t, std::size_t> E320::E320RootDataReader::availableEvents()
    const {
  return {std::get<0>(m_eventMap.front()), std::get<0>(m_eventMap.back()) + 1};
}

ProcessCode E320::E320RootDataReader::read(const AlgorithmContext& ctx) {
  auto it = std::ranges::find_if(m_eventMap, [&](const auto& a) {
    return std::get<0>(a) == ctx.eventNumber;
  });

  if (it == m_eventMap.end()) {
    // Explicitly warn if it happens for the first or last event as that might
    // indicate a human error
    if ((ctx.eventNumber == availableEvents().first) &&
        (ctx.eventNumber == availableEvents().second - 1)) {
      ACTS_WARNING("Reading empty event: " << ctx.eventNumber);
    } else {
      ACTS_DEBUG("Reading empty event: " << ctx.eventNumber);
    }

    m_outputSourceLinks(ctx, {});
    m_outputDetSourceLinksIndices(ctx, {});
    m_outputBpmSourceLinksIndices(ctx, {});

    // Return success flag
    return ProcessCode::SUCCESS;
  }

  // Lock the mutex
  std::lock_guard<std::mutex> lock(m_read_mutex);

  const auto& goInst = *E320::GeometryOptions::instance();

  ACTS_DEBUG("Reading event: " << std::get<0>(*it)
                               << " stored in entries: " << std::get<1>(*it)
                               << " - " << std::get<2>(*it));

  // Create the measurements
  std::vector<Acts::SourceLink> sourceLinks{};
  std::vector<std::size_t> detSourceLinksIndices{};
  std::vector<std::size_t> bpmSourceLinksIndices{};
  EventMetaData eventMetaData{};

  std::size_t eventId = std::get<0>(*it);
  Acts::GeometryIdentifier geoId;
  for (auto entry = std::get<1>(*it); entry < std::get<2>(*it); entry++) {
    m_tree->GetEntry(entry);

    eventMetaData = {m_detEvent->epicsParity, m_detEvent->epicsPID,
                     m_detEvent->epicsDAQNumber};
    for (const auto& staveEv : m_detEvent->st_ev_buffer) {
      for (const auto& chipEv : staveEv.ch_ev_buffer) {
        int sensitiveId = m_geoIdMap.at(chipEv.chip_id);
        if (sensitiveId < m_cfg.minGeoId || sensitiveId > m_cfg.maxGeoId) {
          continue;
        }

        // Apply the Geometry ID convention
        geoId.setSensitive(sensitiveId);

        for (const auto& [hitX, hitY, sizeX, sizeY, size] : chipEv.hits) {
          Acts::Vector2 hitLoc{
              (hitX + 0.5) * goInst.pixelHalfX * 2 - goInst.chipHalfX,
              -(hitY + 0.5) * goInst.pixelHalfY * 2 + goInst.chipHalfY};

          Acts::Vector3 hitGlob = m_cfg.surfaceMap.at(geoId)->localToGlobal(
              ctx.geoContext, hitLoc, Acts::Vector3::UnitX());

          // Estimate error from the cluster size
          Acts::Vector2 stdDev(2 * goInst.pixelHalfX / std::sqrt(12 * size),
                               2 * goInst.pixelHalfY / std::sqrt(12 * size));
          Acts::SquareMatrix2 cov = stdDev.cwiseProduct(stdDev).asDiagonal();

          // Fill the measurement
          SimpleSourceLink ssl(hitLoc, hitGlob, cov, geoId, eventId,
                               sourceLinks.size());
          detSourceLinksIndices.push_back(sourceLinks.size());
          sourceLinks.push_back(Acts::SourceLink(ssl));
        }
      }
    }
    for (const auto& bpmEv : m_detEvent->bpm_ev_buffer) {
      int sensitiveId = m_geoIdMap.at(bpmEv.id);
      if (sensitiveId < m_cfg.minGeoId || sensitiveId > m_cfg.maxGeoId) {
        continue;
      }
      geoId.setSensitive(sensitiveId);

      double quadCenterX = 0;
      double quadCenterY = 0;
      switch (bpmEv.id) {
        case 3156:
          quadCenterX = goInst.quad0CenterShort;
          quadCenterY = goInst.quad0CenterLong;
          break;
        case 3218:
          quadCenterX = goInst.quad1CenterShort;
          quadCenterY = goInst.quad1CenterLong;
          break;
        case 3265:
          quadCenterX = goInst.quad2CenterShort;
          quadCenterY = goInst.quad2CenterLong;
          break;
        case 3315:
          quadCenterX = goInst.quad3CenterShort;
          quadCenterY = goInst.quad3CenterLong;
          break;
        default:
          throw std::runtime_error("Unknown BPM ID");
      }

      Acts::Vector2 hitLoc{bpmEv.meanY + quadCenterY,
                           bpmEv.meanX - quadCenterX};
      Acts::Vector3 hitGlob = m_cfg.surfaceMap.at(geoId)->localToGlobal(
          ctx.geoContext, hitLoc, Acts::Vector3::UnitX());

      // Estimate error
      Acts::Vector2 stdDev(1_mm, 1_mm);
      Acts::SquareMatrix2 cov = stdDev.cwiseProduct(stdDev).asDiagonal();

      // Fill the measurement
      SimpleSourceLink ssl(hitLoc, hitGlob, cov, geoId, eventId,
                           sourceLinks.size());
      bpmSourceLinksIndices.push_back(sourceLinks.size());
      sourceLinks.push_back(Acts::SourceLink(ssl));
    }
  }

  ACTS_DEBUG("Sending " << sourceLinks.size() << " source links");
  ACTS_DEBUG("Sending " << detSourceLinksIndices.size()
                        << " det source link indices");
  ACTS_DEBUG("Sending " << bpmSourceLinksIndices.size()
                        << " bpm source link indices");
  m_outputSourceLinks(ctx, std::move(sourceLinks));
  m_outputDetSourceLinksIndices(ctx, std::move(detSourceLinksIndices));
  m_outputBpmSourceLinksIndices(ctx, std::move(bpmSourceLinksIndices));
  m_outputEventMetaData(ctx, std::move(eventMetaData));

  // Return success flag
  return ProcessCode::SUCCESS;
}
