#include "TrackingPipeline/Io/E320RootDataReader.hpp"

#include <cstddef>

#include <RtypesCore.h>

#include "EpicsFrame.h"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

void E320Io::E320EventFilter::appendAcceleratorState(
    const EpicsFrame& epicsFrame) {
  //--------------------------------------------

  // Energy_BPM_2445_X
  m_state.energy_bpm_2445_x_mean += epicsFrame.energy_bpm_2445_x;
  m_state.energy_bpm_2445_x_stddev +=
      epicsFrame.energy_bpm_2445_x * epicsFrame.energy_bpm_2445_x;
  // Energy_BPM_2445_Y
  m_state.energy_bpm_2445_y_mean += epicsFrame.energy_bpm_2445_y;
  m_state.energy_bpm_2445_y_stddev +=
      epicsFrame.energy_bpm_2445_y * epicsFrame.energy_bpm_2445_y;

  //--------------------------------------------

  // BPM_PB_3156_X
  m_state.bpm_pb_3156_x_mean += epicsFrame.bpm_pb_3156_x;
  m_state.bpm_pb_3156_x_stddev +=
      epicsFrame.bpm_pb_3156_x * epicsFrame.bpm_pb_3156_x;
  // BPM_PB_3156_Y
  m_state.bpm_pb_3156_y_mean += epicsFrame.bpm_pb_3156_y;
  m_state.bpm_pb_3156_y_stddev +=
      epicsFrame.bpm_pb_3156_y * epicsFrame.bpm_pb_3156_y;

  //--------------------------------------------

  // BPM_Q0_3218_X
  m_state.bpm_quad0_3218_x_mean += epicsFrame.bpm_quad0_3218_x;
  m_state.bpm_quad0_3218_x_stddev +=
      epicsFrame.bpm_quad0_3218_x * epicsFrame.bpm_quad0_3218_x;
  // BPM_Q0_3218_Y
  m_state.bpm_quad0_3218_y_mean += epicsFrame.bpm_quad0_3218_y;
  m_state.bpm_quad0_3218_y_stddev +=
      epicsFrame.bpm_quad0_3218_y * epicsFrame.bpm_quad0_3218_y;

  //--------------------------------------------

  // BPM_Q1_3265_X
  m_state.bpm_quad1_3265_x_mean += epicsFrame.bpm_quad1_3265_x;
  m_state.bpm_quad1_3265_x_stddev +=
      epicsFrame.bpm_quad1_3265_x * epicsFrame.bpm_quad1_3265_x;
  // BPM_Q1_3265_Y
  m_state.bpm_quad1_3265_y_mean += epicsFrame.bpm_quad1_3265_y;
  m_state.bpm_quad1_3265_y_stddev +=
      epicsFrame.bpm_quad1_3265_y * epicsFrame.bpm_quad1_3265_y;

  //--------------------------------------------

  // BPM_Q2_3315_X
  m_state.bpm_quad2_3315_x_mean += epicsFrame.bpm_quad2_3315_x;
  m_state.bpm_quad2_3315_x_stddev +=
      epicsFrame.bpm_quad2_3315_x * epicsFrame.bpm_quad2_3315_x;
  // BPM_Q2_3315_Y
  m_state.bpm_quad2_3315_y_mean += epicsFrame.bpm_quad2_3315_y;
  m_state.bpm_quad2_3315_y_stddev +=
      epicsFrame.bpm_quad2_3315_y * epicsFrame.bpm_quad2_3315_y;

  //--------------------------------------------

  // PMT_S20_3060
  m_state.pmt_s20_3060_mean += epicsFrame.pmt_s20_3060;
  m_state.pmt_s20_3060_stddev +=
      epicsFrame.pmt_s20_3060 * epicsFrame.pmt_s20_3060;
  // PMT_S20_3070
  m_state.pmt_s20_3070_mean += epicsFrame.pmt_s20_3070;
  m_state.pmt_s20_3070_stddev +=
      epicsFrame.pmt_s20_3070 * epicsFrame.pmt_s20_3070;
  // PMT_S20_3179
  m_state.pmt_s20_3179_mean += epicsFrame.pmt_s20_3179;
  m_state.pmt_s20_3179_stddev +=
      epicsFrame.pmt_s20_3179 * epicsFrame.pmt_s20_3179;
  // PMT_S20_3350
  m_state.pmt_s20_3350_mean += epicsFrame.pmt_s20_3350;
  m_state.pmt_s20_3350_stddev +=
      epicsFrame.pmt_s20_3350 * epicsFrame.pmt_s20_3350;
  // PMT_S20_3360
  m_state.pmt_s20_3360_mean += epicsFrame.pmt_s20_3360;
  m_state.pmt_s20_3360_stddev +=
      epicsFrame.pmt_s20_3360 * epicsFrame.pmt_s20_3360;

  //--------------------------------------------

  // RADM:LI20:1:CH01:MEAS
  m_state.radm_li20_1_ch01_meas_mean += epicsFrame.radm_li20_1_ch01_meas;
  m_state.radm_li20_1_ch01_meas_stddev +=
      epicsFrame.radm_li20_1_ch01_meas * epicsFrame.radm_li20_1_ch01_meas;
  // XPS:LI20:MC05:M1.RBV
  m_state.xps_li20_mc05_m1_rbv_mean += epicsFrame.xps_li20_mc05_m1_rbv;
  m_state.xps_li20_mc05_m1_rbv_stddev +=
      epicsFrame.xps_li20_mc05_m1_rbv * epicsFrame.xps_li20_mc05_m1_rbv;
  // XPS:LI20:MC05:M2.RBV
  m_state.xps_li20_mc05_m2_rbv_mean += epicsFrame.xps_li20_mc05_m2_rbv;
  m_state.xps_li20_mc05_m2_rbv_stddev +=
      epicsFrame.xps_li20_mc05_m2_rbv * epicsFrame.xps_li20_mc05_m2_rbv;

  m_stateSize++;
}

void E320Io::E320EventFilter::finalizeAcceleratorState() {
  //--------------------------------------------

  // Energy_BPM_2445_X
  m_state.energy_bpm_2445_x_mean /= m_stateSize;
  m_state.energy_bpm_2445_x_stddev /= m_stateSize;
  m_state.energy_bpm_2445_x_stddev = std::sqrt(
      m_state.energy_bpm_2445_x_stddev -
      m_state.energy_bpm_2445_x_mean * m_state.energy_bpm_2445_x_mean);
  // Energy_BPM_2445_Y
  m_state.energy_bpm_2445_y_mean /= m_stateSize;
  m_state.energy_bpm_2445_y_stddev /= m_stateSize;
  m_state.energy_bpm_2445_y_stddev = std::sqrt(
      m_state.energy_bpm_2445_y_stddev -
      m_state.energy_bpm_2445_y_mean * m_state.energy_bpm_2445_y_mean);

  //--------------------------------------------

  // BPM_PB_3156_X
  m_state.bpm_pb_3156_x_mean /= m_stateSize;
  m_state.bpm_pb_3156_x_stddev /= m_stateSize;
  m_state.bpm_pb_3156_x_stddev =
      std::sqrt(m_state.bpm_pb_3156_x_stddev -
                m_state.bpm_pb_3156_x_mean * m_state.bpm_pb_3156_x_mean);
  // BPM_PB_3156_Y
  m_state.bpm_pb_3156_y_mean /= m_stateSize;
  m_state.bpm_pb_3156_y_stddev /= m_stateSize;
  m_state.bpm_pb_3156_y_stddev =
      std::sqrt(m_state.bpm_pb_3156_y_stddev -
                m_state.bpm_pb_3156_y_mean * m_state.bpm_pb_3156_y_mean);

  //--------------------------------------------

  // BPM_Q0_3218_X
  m_state.bpm_quad0_3218_x_mean /= m_stateSize;
  m_state.bpm_quad0_3218_x_stddev /= m_stateSize;
  m_state.bpm_quad0_3218_x_stddev =
      std::sqrt(m_state.bpm_quad0_3218_x_stddev -
                m_state.bpm_quad0_3218_x_mean * m_state.bpm_quad0_3218_x_mean);
  // BPM_Q0_3218_Y
  m_state.bpm_quad0_3218_y_mean /= m_stateSize;
  m_state.bpm_quad0_3218_y_stddev /= m_stateSize;
  m_state.bpm_quad0_3218_y_stddev =
      std::sqrt(m_state.bpm_quad0_3218_y_stddev -
                m_state.bpm_quad0_3218_y_mean * m_state.bpm_quad0_3218_y_mean);

  //--------------------------------------------

  // BPM_Q1_3265_X
  m_state.bpm_quad1_3265_x_mean /= m_stateSize;
  m_state.bpm_quad1_3265_x_stddev /= m_stateSize;
  m_state.bpm_quad1_3265_x_stddev =
      std::sqrt(m_state.bpm_quad1_3265_x_stddev -
                m_state.bpm_quad1_3265_x_mean * m_state.bpm_quad1_3265_x_mean);
  // BPM_Q1_3265_Y
  m_state.bpm_quad1_3265_y_mean /= m_stateSize;
  m_state.bpm_quad1_3265_y_stddev /= m_stateSize;
  m_state.bpm_quad1_3265_y_stddev =
      std::sqrt(m_state.bpm_quad1_3265_y_stddev -
                m_state.bpm_quad1_3265_y_mean * m_state.bpm_quad1_3265_y_mean);

  //--------------------------------------------

  // BPM_Q2_3315_X
  m_state.bpm_quad2_3315_x_mean /= m_stateSize;
  m_state.bpm_quad2_3315_x_stddev /= m_stateSize;
  m_state.bpm_quad2_3315_x_stddev =
      std::sqrt(m_state.bpm_quad2_3315_x_stddev -
                m_state.bpm_quad2_3315_x_mean * m_state.bpm_quad2_3315_x_mean);
  // BPM_Q2_3315_Y
  m_state.bpm_quad2_3315_y_mean /= m_stateSize;
  m_state.bpm_quad2_3315_y_stddev /= m_stateSize;
  m_state.bpm_quad2_3315_y_stddev =
      std::sqrt(m_state.bpm_quad2_3315_y_stddev -
                m_state.bpm_quad2_3315_y_mean * m_state.bpm_quad2_3315_y_mean);

  //--------------------------------------------

  // PMT_S20_3060
  m_state.pmt_s20_3060_mean /= m_stateSize;
  m_state.pmt_s20_3060_stddev /= m_stateSize;
  m_state.pmt_s20_3060_stddev =
      std::sqrt(m_state.pmt_s20_3060_stddev -
                m_state.pmt_s20_3060_mean * m_state.pmt_s20_3060_mean);
  // PMT_S20_3070
  m_state.pmt_s20_3070_mean /= m_stateSize;
  m_state.pmt_s20_3070_stddev /= m_stateSize;
  m_state.pmt_s20_3070_stddev =
      std::sqrt(m_state.pmt_s20_3070_stddev -
                m_state.pmt_s20_3070_mean * m_state.pmt_s20_3070_mean);
  // PMT_S20_3179
  m_state.pmt_s20_3179_mean /= m_stateSize;
  m_state.pmt_s20_3179_stddev /= m_stateSize;
  m_state.pmt_s20_3179_stddev =
      std::sqrt(m_state.pmt_s20_3179_stddev -
                m_state.pmt_s20_3179_mean * m_state.pmt_s20_3179_mean);
  // PMT_S20_3350
  m_state.pmt_s20_3350_mean /= m_stateSize;
  m_state.pmt_s20_3350_stddev /= m_stateSize;
  m_state.pmt_s20_3350_stddev =
      std::sqrt(m_state.pmt_s20_3350_stddev -
                m_state.pmt_s20_3350_mean * m_state.pmt_s20_3350_mean);
  // PMT_S20_3360
  m_state.pmt_s20_3360_mean /= m_stateSize;
  m_state.pmt_s20_3360_stddev /= m_stateSize;
  m_state.pmt_s20_3360_stddev =
      std::sqrt(m_state.pmt_s20_3360_stddev -
                m_state.pmt_s20_3360_mean * m_state.pmt_s20_3360_mean);

  //--------------------------------------------

  // RADM:LI20:1:CH01:MEAS
  m_state.radm_li20_1_ch01_meas_mean /= m_stateSize;
  m_state.radm_li20_1_ch01_meas_stddev /= m_stateSize;
  m_state.radm_li20_1_ch01_meas_stddev = std::sqrt(
      m_state.radm_li20_1_ch01_meas_stddev -
      m_state.radm_li20_1_ch01_meas_mean * m_state.radm_li20_1_ch01_meas_mean);
  // XPS:LI20:MC05:M1.RBV
  m_state.xps_li20_mc05_m1_rbv_mean /= m_stateSize;
  m_state.xps_li20_mc05_m1_rbv_stddev /= m_stateSize;
  m_state.xps_li20_mc05_m1_rbv_stddev = std::sqrt(
      m_state.xps_li20_mc05_m1_rbv_stddev -
      m_state.xps_li20_mc05_m1_rbv_mean * m_state.xps_li20_mc05_m1_rbv_mean);
  // XPS:LI20:MC05:M2.RBV
  m_state.xps_li20_mc05_m2_rbv_mean /= m_stateSize;
  m_state.xps_li20_mc05_m2_rbv_stddev /= m_stateSize;
  m_state.xps_li20_mc05_m2_rbv_stddev = std::sqrt(
      m_state.xps_li20_mc05_m2_rbv_stddev -
      m_state.xps_li20_mc05_m2_rbv_mean * m_state.xps_li20_mc05_m2_rbv_mean);
}

bool E320Io::E320EventFilter::checkCuts(const DetectorEvent& detEvent) {
  // Discriminate events based on the turnaround time
  if ((detEvent.ts_end - detEvent.ts_begin) / 1e6 > m_cfg.maxTurnaroundTime) {
    return false;
  }

  for (const auto& staveEvent : detEvent.st_ev_buffer) {
    for (const auto& chipEvent : staveEvent.ch_ev_buffer) {
      if (chipEvent.is_busy_violation || chipEvent.is_flushed_incomplete ||
          chipEvent.is_strobe_extended || chipEvent.is_busy_transition) {
        return false;
      }
      if (chipEvent.end_of_run || chipEvent.overflow || chipEvent.timeout ||
          chipEvent.header_error || chipEvent.decoder_10b8b_error ||
          chipEvent.event_oversize_error) {
        return false;
      }
    }
  }
  //--------------------------------------------

  // Energy_BPM_2445_X
  if (std::abs(detEvent.epics_frame.energy_bpm_2445_x -
               m_state.energy_bpm_2445_x_mean) /
          m_state.energy_bpm_2445_x_stddev >
      m_cfg.nStdDevs) {
    return false;
  }
  // Energy_BPM_2445_Y
  if (std::abs(detEvent.epics_frame.energy_bpm_2445_y -
               m_state.energy_bpm_2445_y_mean) /
          m_state.energy_bpm_2445_y_stddev >
      m_cfg.nStdDevs) {
    return false;
  }

  //--------------------------------------------

  // BPM_PB_3156_X
  if (std::abs(detEvent.epics_frame.bpm_pb_3156_x -
               m_state.bpm_pb_3156_x_mean) /
          m_state.bpm_pb_3156_x_stddev >
      m_cfg.nStdDevs) {
    return false;
  }
  // BPM_PB_3156_Y
  if (std::abs(detEvent.epics_frame.bpm_pb_3156_y -
               m_state.bpm_pb_3156_y_mean) /
          m_state.bpm_pb_3156_y_stddev >
      m_cfg.nStdDevs) {
    return false;
  }

  //--------------------------------------------

  // BPM_Q0_3218_X
  if (std::abs(detEvent.epics_frame.bpm_quad0_3218_x -
               m_state.bpm_quad0_3218_x_mean) /
          m_state.bpm_quad0_3218_x_stddev >
      m_cfg.nStdDevs) {
    return false;
  }
  // BPM_Q0_3218_Y
  if (std::abs(detEvent.epics_frame.bpm_quad0_3218_y -
               m_state.bpm_quad0_3218_y_mean) /
          m_state.bpm_quad0_3218_y_stddev >
      m_cfg.nStdDevs) {
    return false;
  }

  //--------------------------------------------

  // BPM_Q1_3265_X
  if (std::abs(detEvent.epics_frame.bpm_quad1_3265_x -
               m_state.bpm_quad1_3265_x_mean) /
          m_state.bpm_quad1_3265_x_stddev >
      m_cfg.nStdDevs) {
    return false;
  }
  // BPM_Q1_3265_Y
  if (std::abs(detEvent.epics_frame.bpm_quad1_3265_y -
               m_state.bpm_quad1_3265_y_mean) /
          m_state.bpm_quad1_3265_y_stddev >
      m_cfg.nStdDevs) {
    return false;
  }

  //--------------------------------------------

  // BPM_Q2_3315_X
  if (std::abs(detEvent.epics_frame.bpm_quad2_3315_x -
               m_state.bpm_quad2_3315_x_mean) /
          m_state.bpm_quad2_3315_x_stddev >
      m_cfg.nStdDevs) {
    return false;
  }
  // BPM_Q2_3315_Y
  if (std::abs(detEvent.epics_frame.bpm_quad2_3315_y -
               m_state.bpm_quad2_3315_y_mean) /
          m_state.bpm_quad2_3315_y_stddev >
      m_cfg.nStdDevs) {
    return false;
  }

  //--------------------------------------------

  // PMT_S20_3060
  if (std::abs(detEvent.epics_frame.pmt_s20_3060 - m_state.pmt_s20_3060_mean) /
          m_state.pmt_s20_3060_stddev >
      m_cfg.nStdDevs) {
    return false;
  }
  // PMT_S20_3070
  if (std::abs(detEvent.epics_frame.pmt_s20_3070 - m_state.pmt_s20_3070_mean) /
          m_state.pmt_s20_3070_stddev >
      m_cfg.nStdDevs) {
    return false;
  }
  // PMT_S20_3179
  if (std::abs(detEvent.epics_frame.pmt_s20_3179 - m_state.pmt_s20_3179_mean) /
          m_state.pmt_s20_3179_stddev >
      m_cfg.nStdDevs) {
    return false;
  }
  // PMT_S20_3350
  if (std::abs(detEvent.epics_frame.pmt_s20_3350 - m_state.pmt_s20_3350_mean) /
          m_state.pmt_s20_3350_stddev >
      m_cfg.nStdDevs) {
    return false;
  }
  // PMT_S20_3360
  if (std::abs(detEvent.epics_frame.pmt_s20_3360 - m_state.pmt_s20_3360_mean) /
          m_state.pmt_s20_3360_stddev >
      m_cfg.nStdDevs) {
    return false;
  }

  //--------------------------------------------

  // RADM:LI20:1:CH01:MEAS
  if (std::abs(detEvent.epics_frame.radm_li20_1_ch01_meas -
               m_state.radm_li20_1_ch01_meas_mean) /
          m_state.radm_li20_1_ch01_meas_stddev >
      m_cfg.nStdDevs) {
    return false;
  }
  // XPS:LI20:MC05:M1.RBV
  if (std::abs(detEvent.epics_frame.xps_li20_mc05_m1_rbv -
               m_state.xps_li20_mc05_m1_rbv_mean) /
          m_state.xps_li20_mc05_m1_rbv_stddev >
      m_cfg.nStdDevs) {
    return false;
  }
  // XPS:LI20:MC05:M2.RBV
  if (std::abs(detEvent.epics_frame.xps_li20_mc05_m2_rbv -
               m_state.xps_li20_mc05_m2_rbv_mean) /
          m_state.xps_li20_mc05_m2_rbv_stddev >
      m_cfg.nStdDevs) {
    return false;
  }

  return true;
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
  m_chain->SetBranchAddress("eventId", &m_eventId);
  m_chain->SetBranchAddress(m_cfg.eventKey.c_str(), &m_detEvent);

  // Add the files to the chain
  for (const auto& path : m_cfg.filePaths) {
    m_chain->Add(path.c_str());
  }
  auto nEntries = static_cast<std::size_t>(m_chain->GetEntries());

  // Add the first entry
  m_chain->GetEntry(m_cfg.skip);
  m_eventMap.emplace_back(m_eventId, 0, 0);
  if (m_cfg.dataFilter != nullptr) {
    m_cfg.dataFilter->appendAcceleratorState(m_detEvent->epics_frame);
  }

  // Go through all entries and store the position of the events
  for (std::size_t i = m_cfg.skip + 1; i < nEntries; ++i) {
    m_chain->GetEntry(i);
    if (m_cfg.dataFilter != nullptr) {
      m_cfg.dataFilter->appendAcceleratorState(m_detEvent->epics_frame);
    }
    if (m_eventId != std::get<0>(m_eventMap.back())) {
      std::get<2>(m_eventMap.back()) = i;
      m_eventMap.emplace_back(m_eventId, i, i);
    }
  }

  if (m_cfg.dataFilter != nullptr) {
    m_cfg.dataFilter->finalizeAcceleratorState();
  }

  // Sort by event id
  std::ranges::sort(m_eventMap, [](const auto& a, const auto& b) {
    return std::get<0>(a) < std::get<0>(b);
  });

  std::get<2>(m_eventMap.back()) = nEntries;

  // Add the data branch
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

    // Check if the event is good
    if (m_cfg.dataFilter != nullptr &&
        !m_cfg.dataFilter->checkCuts(*m_detEvent)) {
      return ProcessCode::SUCCESS;
    }

    for (const auto& staveEv : m_detEvent->st_ev_buffer) {
      for (const auto& chipEv : staveEv.ch_ev_buffer) {
        // Apply the Geometry ID convention
        std::size_t geoIdVal = static_cast<std::size_t>(chipEv.chip_id + 1);
        Acts::GeometryIdentifier geoId;
        geoId.setSensitive(geoIdVal);

        for (const auto& [hitX, hitY, sizeX, sizeY, size] : chipEv.hits) {
          if (m_cfg.clusterFilter != nullptr &&
              !m_cfg.clusterFilter->operator()(hitX, hitY)) {
            continue;
          }
          if (size > 4) {
            continue;
          }

          Acts::Vector2 hitLoc{hitX, hitY};

          // TODO: Estimate from the simulation
          // Estimate error from the cluster size
          double errX = m_gOpt.pixelSizeY * sizeY;
          double errY = m_gOpt.pixelSizeX * sizeX;
          /*double errX = m_gOpt.pixelSizeY / std::sqrt(12 * size);*/
          /*double errY = m_gOpt.pixelSizeX / std::sqrt(12 * size);*/

          Acts::Vector2 stdDev(errX, errY);
          Acts::SquareMatrix2 cov = stdDev.cwiseProduct(stdDev).asDiagonal();

          // Fill the measurement
          SimpleSourceLink ssl(hitLoc, cov, geoId, eventId, sourceLinks.size());
          sourceLinks.push_back(Acts::SourceLink(ssl));
        }
      }
    }
  }

  m_outputSourceLinks(context, std::move(sourceLinks));

  // Return success flag
  return ProcessCode::SUCCESS;
}
