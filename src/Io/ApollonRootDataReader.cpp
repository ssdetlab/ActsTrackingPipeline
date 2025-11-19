#include "TrackingPipeline/Io/ApollonRootDataReader.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

namespace ag = ApollonGeometry;

ApollonIo::ApollonRootDataReader::ApollonRootDataReader(
    const Config& config, Acts::Logging::Level level)
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

std::pair<std::size_t, std::size_t>
ApollonIo::ApollonRootDataReader::availableEvents() const {
  return {std::get<0>(m_eventMap.front()), std::get<0>(m_eventMap.back()) + 1};
}

ProcessCode ApollonIo::ApollonRootDataReader::read(
    const AlgorithmContext& context) {
  auto it = std::ranges::find_if(m_eventMap, [&](const auto& a) {
    return std::get<0>(a) == context.eventNumber;
  });

  if (it == m_eventMap.end()) {
    // explicitly warn if it happens for the first or last event as that might
    // indicate a human error
    if ((context.eventNumber == availableEvents().first) &&
        (context.eventNumber == availableEvents().second - 1)) {
      ACTS_WARNING("Reading empty event: " << context.eventNumber);
    } else {
      ACTS_DEBUG("Reading empty event: " << context.eventNumber);
    }

    m_outputSourceLinks(context, {});

    // Return success flag
    return ProcessCode::SUCCESS;
  }

  // lock the mutex
  std::lock_guard<std::mutex> lock(m_read_mutex);

  const auto& goInst = *ag::GeometryOptions::instance();

  ACTS_DEBUG("Reading event: " << std::get<0>(*it)
                               << " stored in entries: " << std::get<1>(*it)
                               << " - " << std::get<2>(*it));

  // Create the measurements
  std::vector<Acts::SourceLink> sourceLinks{};
  std::size_t eventId = std::get<0>(*it);
  Acts::GeometryIdentifier geoId;
  for (auto entry = std::get<1>(*it); entry < std::get<2>(*it); entry++) {
    m_chain->GetEntry(entry);

    for (const auto& staveEv : m_detEvent->st_ev_buffer) {
      for (const auto& chipEv : staveEv.ch_ev_buffer) {
        // Apply the Geometry ID convention
        std::size_t geoIdVal =
            (staveEv.stave_id == 1)
                ? static_cast<std::size_t>(10 + chipEv.chip_id)
                : static_cast<std::size_t>(20 + chipEv.chip_id);
        if (geoIdVal < m_cfg.minGeoId || geoIdVal > m_cfg.maxGeoId) {
          continue;
        }
        geoId.setSensitive(geoIdVal);

        for (const auto& [hitX, hitY, sizeX, sizeY, size] : chipEv.hits) {
          Acts::Vector2 hitLoc{
              (hitX + 0.5) * goInst.pixelHalfX * 2 - goInst.chipHalfX,
              (hitY + 0.5) * goInst.pixelHalfY * 2 - goInst.chipHalfY};
          Acts::Vector3 hitGlob = m_cfg.surfaceMap.at(geoId)->localToGlobal(
              context.geoContext, hitLoc, Acts::Vector3::UnitX());

          // Estimate error from the cluster size
          Acts::Vector2 stdDev(2 * goInst.pixelHalfX / std::sqrt(12 * size),
                               2 * goInst.pixelHalfY / std::sqrt(12 * size));
          Acts::SquareMatrix2 cov = stdDev.cwiseProduct(stdDev).asDiagonal();

          // Fill the measurement
          SimpleSourceLink ssl(hitLoc, hitGlob, cov, geoId, eventId,
                               sourceLinks.size());
          sourceLinks.push_back(Acts::SourceLink(ssl));
        }
      }
    }
  }

  ACTS_DEBUG("Sending " << sourceLinks.size() << " measurements");
  m_outputSourceLinks(context, std::move(sourceLinks));

  // Return success flag
  return ProcessCode::SUCCESS;
}
