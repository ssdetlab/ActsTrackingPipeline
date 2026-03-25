#include "TrackingPipeline/Io/E320RootDataReader.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <stdexcept>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

E320::E320RootDataReader::E320RootDataReader(const Config& config,
                                             Acts::Logging::Level level)
    : IReader(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  // m_chain = new TChain(m_cfg.treeName.c_str());
  m_file = new TFile(m_cfg.filePaths.at(0).c_str());
  m_chain = m_file->Get<TTree>(m_cfg.treeName.c_str());

  if (m_cfg.filePaths.empty()) {
    throw std::invalid_argument("Missing input filenames");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);

  // Set the branches
  m_chain->SetBranchAddress("eventId", &m_eventId);
  m_chain->SetBranchAddress(m_cfg.eventKey.c_str(), &m_detEvent);

  // Add the files to the chain
  // for (const auto& path : m_cfg.filePaths) {
  //   m_chain->Add(path.c_str());
  // }
  m_chain->SetBranchStatus("*", false);
  if (!m_chain->GetBranch("eventId")) {
    throw std::invalid_argument("Missing eventId branch");
  }
  m_chain->SetBranchStatus("eventId", true);
  auto nEntries = static_cast<std::size_t>(m_chain->GetEntries());

  // Add the first entry
  m_chain->GetEntry(0);
  m_eventMap.emplace_back(m_eventId, 0, 0);

  // Go through all entries and store the position of the events
  for (std::size_t i = 1; i < nEntries; ++i) {
    m_chain->GetEntry(i);
    if (m_eventId != std::get<0>(m_eventMap.back())) {
      std::get<2>(m_eventMap.back()) = i;
      m_eventMap.emplace_back(m_eventId, i, i);
    }
  }

  // Sort by event id
  std::ranges::sort(m_eventMap, [](const auto& a, const auto& b) {
    return std::get<0>(a) < std::get<0>(b);
  });

  std::get<2>(m_eventMap.back()) = nEntries;

  // Add the data branch
  m_chain->SetBranchStatus("*", true);

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
    // explicitly warn if it happens for the first or last event as that might
    // indicate a human error
    if ((ctx.eventNumber == availableEvents().first) &&
        (ctx.eventNumber == availableEvents().second - 1)) {
      ACTS_WARNING("Reading empty event: " << ctx.eventNumber);
    } else {
      ACTS_DEBUG("Reading empty event: " << ctx.eventNumber);
    }

    m_outputSourceLinks(ctx, {});

    // Return success flag
    return ProcessCode::SUCCESS;
  }

  // lock the mutex
  std::lock_guard<std::mutex> lock(m_read_mutex);

  const auto& goInst = *E320::GeometryOptions::instance();

  ACTS_DEBUG("Reading event: " << std::get<0>(*it)
                               << " stored in entries: " << std::get<1>(*it)
                               << " - " << std::get<2>(*it));

  std::unordered_map<int, std::pair<double, double>> clStdDevs = {
      {1, {0.00395262, 0.00333208}},
      {2, {0.00444654, 0.00473498}},
      {3, {0.00417356, 0.00431227}},
      {4, {0.00509384, 0.00490848}}};

  // Create the measurements
  std::vector<Acts::SourceLink> sourceLinks{};
  std::size_t eventId = std::get<0>(*it);
  Acts::GeometryIdentifier geoId;
  for (auto entry = std::get<1>(*it); entry < std::get<2>(*it); entry++) {
    m_chain->GetEntry(entry);

    for (const auto& staveEv : m_detEvent->st_ev_buffer) {
      for (const auto& chipEv : staveEv.ch_ev_buffer) {
        int sensitiveId = m_geoIdMap.at(chipEv.chip_id);
        if (sensitiveId < m_cfg.minGeoId || sensitiveId > m_cfg.maxGeoId) {
          continue;
        }

        // Apply the Geometry ID convention
        geoId.setSensitive(sensitiveId);

        if (chipEv.hits.size() > 1000) {
          continue;
        }
        for (const auto& [hitX, hitY, sizeX, sizeY, size] : chipEv.hits) {
          Acts::Vector2 hitLoc{
              (hitX + 0.5) * goInst.pixelHalfX * 2 - goInst.chipHalfX,
              -(hitY + 0.5) * goInst.pixelHalfY * 2 + goInst.chipHalfY};

          Acts::Vector3 hitGlob = m_cfg.surfaceMap.at(geoId)->localToGlobal(
              ctx.geoContext, hitLoc, Acts::Vector3::UnitX());

          // Estimate error from the cluster size
          Acts::Vector2 stdDev(clStdDevs.at(size).first,
                               clStdDevs.at(size).second);
          // Acts::Vector2 stdDev(2 * goInst.pixelHalfX / std::sqrt(12 * size),
          //                      2 * goInst.pixelHalfY / std::sqrt(12 * size));
          Acts::SquareMatrix2 cov = stdDev.cwiseProduct(stdDev).asDiagonal();

          // Fill the measurement
          SimpleSourceLink ssl(hitLoc, hitGlob, cov, geoId, eventId,
                               sourceLinks.size());

          Acts::Vector2 stdDevInf(2 * goInst.pixelHalfX * sizeX,
                                  2 * goInst.pixelHalfY * sizeY);
          Acts::SquareMatrix2 covInf =
              stdDevInf.cwiseProduct(stdDevInf).asDiagonal();
          ssl.setCovarianceInf(covInf);

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

      double quadCenterX;
      double quadCenterY;
      switch (bpmEv.id) {
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

      Acts::Vector2 hitLoc{bpmEv.meanX + quadCenterX,
                           -(bpmEv.meanY + quadCenterY)};
      Acts::Vector3 hitGlob = m_cfg.surfaceMap.at(geoId)->localToGlobal(
          ctx.geoContext, hitLoc, Acts::Vector3::UnitX());

      // Estimate error from the cluster size
      Acts::Vector2 stdDev(bpmEv.stdDevX, bpmEv.stdDevY);
      Acts::SquareMatrix2 cov = stdDev.cwiseProduct(stdDev).asDiagonal();

      // Fill the measurement
      SimpleSourceLink ssl(hitLoc, hitGlob, cov, geoId, eventId,
                           sourceLinks.size());
      Acts::SquareMatrix2 covInf = cov;
      ssl.setCovarianceInf(covInf);

      sourceLinks.push_back(Acts::SourceLink(ssl));
    }
  }

  ACTS_DEBUG("Sending " << sourceLinks.size() << " measurements");
  m_outputSourceLinks(ctx, std::move(sourceLinks));

  // Return success flag
  return ProcessCode::SUCCESS;
}
