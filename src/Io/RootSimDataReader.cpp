#include "TrackingPipeline/Io/RootSimDataReader.hpp"

// Set branch addresses for the TChain
template <typename K, typename T>
void setBranches(TChain* chain, const K& keys,
                 std::unordered_map<std::string_view, T>& columns) {
  T value = 0;
  for (auto key : keys) {
    columns.insert({key, value});
  }
  for (auto key : keys) {
    chain->SetBranchAddress(key, &columns.at(key));
  }
};

RootSimDataReader::RootSimDataReader(const Config& config,
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
  m_outputSimClusters.initialize(m_cfg.outputSimClusters);

  // Set the branches
  setBranches(m_chain, m_cfg.vVector3Keys, m_vVector3Columns);
  setBranches(m_chain, m_cfg.vector3Keys, m_vector3Columns);
  setBranches(m_chain, m_cfg.vLorentzKeys, m_vLorentzColumns);
  setBranches(m_chain, m_cfg.vIntKeys, m_vIntColumns);
  setBranches(m_chain, m_cfg.vDoubleKeys, m_vDoubleColumns);
  setBranches(m_chain, m_cfg.intKeys, m_intColumns);

  // Add the files to the chain
  for (const auto& path : m_cfg.filePaths) {
    m_chain->Add(path.c_str());
  }

  // Disable all branches and only enable event-id for a first scan of the
  // file
  m_chain->SetBranchStatus("*", false);
  if (!m_chain->GetBranch("eventId")) {
    throw std::invalid_argument("Missing eventId branch");
  }
  m_chain->SetBranchStatus("eventId", true);

  auto nEntries = static_cast<std::size_t>(m_chain->GetEntries());

  // Add the first entry
  m_chain->GetEntry(0);
  m_eventMap.emplace_back(m_intColumns.at("eventId"), 0, 0);

  // Go through all entries and store the position of the events
  for (std::size_t i = 1; i < nEntries; ++i) {
    m_chain->GetEntry(i);
    const auto evtId = m_intColumns.at("eventId");

    if (evtId != std::get<0>(m_eventMap.back())) {
      std::get<2>(m_eventMap.back()) = i;
      m_eventMap.emplace_back(evtId, i, i);
    }
  }
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

std::pair<std::size_t, std::size_t> RootSimDataReader::availableEvents() const {
  return {std::get<0>(m_eventMap.front()), std::get<0>(m_eventMap.back()) + 1};
}

ProcessCode RootSimDataReader::read(const AlgorithmContext& context) {
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
    m_outputSimClusters(context, {});

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
  SimClusters clusters{};
  for (auto entry = std::get<1>(*it); entry < std::get<2>(*it); entry++) {
    m_chain->GetEntry(entry);
    prepareMeasurements(context, &sourceLinks, &clusters);
  }

  m_outputSourceLinks(context, std::move(sourceLinks));
  m_outputSimClusters(context, std::move(clusters));

  // Return success flag
  return ProcessCode::SUCCESS;
}
