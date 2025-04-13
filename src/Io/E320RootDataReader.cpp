#include "TrackingPipeline/Io/E320RootDataReader.hpp"

#include <cstddef>

#include <RtypesCore.h>

void E320Io::E320RootDataReader::bfsClustering(std::vector<Hit>& clusterPixels,
                                               Hit& pivot,
                                               std::set<Hit>& pixels) {
  std::deque<std::pair<int, int>> queue{pivot};
  while (!queue.empty()) {
    auto current = queue.front();
    queue.pop_front();
    auto currentPtr = pixels.find(current);

    if (currentPtr != pixels.end()) {
      // add this pixel to the cluster
      clusterPixels.push_back(current);
      // remove pixel from live pixels list
      pixels.erase(currentPtr);

      std::vector<Hit> neighbours{{current.first + 1, current.second},
                                  {current.first - 1, current.second},
                                  {current.first, current.second + 1},
                                  {current.first, current.second - 1}};
      for (const auto& n : neighbours) {
        auto nptr = pixels.find(n);
        if (nptr != pixels.end()) {
          queue.push_back(*nptr);
        }
      }
    }
  }
}

std::vector<E320Io::E320RootDataReader::Cluster>
E320Io::E320RootDataReader::getClusters(std::set<Hit>& hits) {
  std::vector<Cluster> clusters;

  while (!hits.empty()) {
    auto pixel = *hits.begin();
    std::vector<Hit> clusterPixels;
    bfsClustering(clusterPixels, pixel, hits);

    std::size_t xMin = 1024;
    std::size_t xMax = 0;
    std::size_t yMin = 512;
    std::size_t yMax = 0;
    double xCenter = 0;
    double yCenter = 0;
    for (auto pix : clusterPixels) {
      xCenter += pix.first;
      yCenter += pix.second;

      if (pix.first <= xMin) {
        xMin = pix.first;
      }
      if (pix.first >= xMax) {
        xMax = pix.first;
      }
      if (pix.second <= yMin) {
        yMin = pix.second;
      }
      if (pix.second >= yMax) {
        yMax = pix.second;
      }
    }
    std::size_t dx = xMax - xMin + 1;
    std::size_t dy = yMax - yMin + 1;
    std::size_t cl_size = clusterPixels.size();
    xCenter /= clusterPixels.size();
    yCenter /= clusterPixels.size();

    clusters.push_back(Cluster{xCenter, yCenter, dx, dy, cl_size});
  }

  return clusters;
}

E320Io::E320RootDataReader::E320RootDataReader(const Config& config,
                                               Acts::Logging::Level level)
    : IReader(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  m_chain = new TChain(m_cfg.treeName.c_str());

  if (m_cfg.filePaths.empty()) {
    throw std::invalid_argument("Missing input filenames");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);

  // Set the branches
  m_chain->SetBranchAddress(m_cfg.eventKey.c_str(), &m_detEvent);
  m_chain->SetBranchAddress("eventId", &m_eventId);

  // Add the files to the chain
  for (const auto& path : m_cfg.filePaths) {
    m_chain->Add(path.c_str());
  }

  // Disable all branches and only enable event-id for a first scan of the
  // file
  //  m_chain->SetBranchStatus("*", false);
  //  if (!m_chain->GetBranch("eventId")) {
  //    throw std::invalid_argument("Missing eventId branch");
  //  }
  //  m_chain->SetBranchStatus("eventId", true);

  auto nEntries = static_cast<std::size_t>(m_chain->GetEntries());

  // Add the first entry
  m_chain->GetEntry(0);
  m_eventMap.emplace_back(m_eventId, 0, 0);

  // Go through all entries and store the position of the events
  for (std::size_t i = 1; i < nEntries; ++i) {
    m_chain->GetEntry(i);
    std::cout << "eventID = " << m_eventId << "\n";
    if (m_eventId != std::get<0>(m_eventMap.back())) {
      std::get<2>(m_eventMap.back()) = i;
      m_eventMap.emplace_back(m_eventId, i, i);
    }
  }
  std::cout << "Events " << m_eventMap.size() << "\n";

  // Sort by event id
  std::ranges::sort(m_eventMap, [](const auto& a, const auto& b) {
    return std::get<0>(a) < std::get<0>(b);
  });

  std::get<2>(m_eventMap.back()) = nEntries;

  // Re-Enable all branches
  m_chain->SetBranchStatus("*", true);
  ACTS_DEBUG("Event range: " << availableEvents().first << " - "
                             << availableEvents().second);
}

std::pair<std::size_t, std::size_t>
E320Io::E320RootDataReader::availableEvents() const {
  return {std::get<0>(m_eventMap.front()), std::get<0>(m_eventMap.back()) + 1};
}

ProcessCode E320Io::E320RootDataReader::read(const AlgorithmContext& context) {
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

  ACTS_DEBUG("Reading event: " << std::get<0>(*it)
                               << " stored in entries: " << std::get<1>(*it)
                               << " - " << std::get<2>(*it));

  // Create the measurements
  std::vector<Acts::SourceLink> sourceLinks{};
  std::size_t eventId = std::get<0>(*it);
  for (auto entry = std::get<1>(*it); entry < std::get<2>(*it); entry++) {
    m_chain->GetEntry(entry);

    for (const auto& staveEv : m_detEvent->st_ev_buffer) {
      for (const auto& chipEv : staveEv.ch_ev_buffer) {
        std::set<std::pair<std::size_t, std::size_t>> hits{chipEv.hits.begin(),
                                                           chipEv.hits.end()};
        // Apply the Geometry ID convention
        std::size_t geoIdVal = static_cast<std::size_t>(chipEv.chip_id + 1);
        Acts::GeometryIdentifier geoId;
        geoId.setSensitive(geoIdVal);

        auto clusters = getClusters(hits);
        for (const auto& [hitX, hitY, sizeX, sizeY, size] : clusters) {
          // Convert cluster center to mm
          double xPix = (hitX + 0.5) * m_gOpt.pixelSizeY;
          double yPix = (hitY + 0.5) * m_gOpt.pixelSizeX;

          // PI rotation
          xPix *= -1;
          yPix *= -1;

          // Move the origin to the chip's center
          double xPixLoc = xPix + m_gOpt.chipSizeY / 2;
          double yPixLoc = yPix + m_gOpt.chipSizeX / 2;

          Acts::Vector2 hitLoc{xPixLoc, yPixLoc};

          // Estimate error from the cluster size
          double errX = m_gOpt.pixelSizeY / std::sqrt(12 * size);
          double errY = m_gOpt.pixelSizeX / std::sqrt(12 * size);

          Acts::Vector2 stdDev(errX, errY);
          Acts::SquareMatrix2 cov = stdDev.cwiseProduct(stdDev).asDiagonal();

          // Fill the measurement
          SimpleSourceLink ssl(hitLoc, cov, geoId, eventId, sourceLinks.size());
        }
      }
    }
  }

  m_outputSourceLinks(context, std::move(sourceLinks));

  // Return success flag
  return ProcessCode::SUCCESS;
}
