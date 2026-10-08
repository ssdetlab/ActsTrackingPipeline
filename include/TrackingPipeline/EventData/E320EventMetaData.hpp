#pragma once

#include <cstddef>

namespace E320 {

/// @brief Event meta data struct
struct EventMetaData {
  /// MOSAIC trigger number
  std::size_t eudaqTrgN;
  /// EUDAQ DAQ number
  std::size_t eudaqDAQNumber;
  /// EUDAQ run start timestamp
  std::size_t eudaqRunStartTs;
  /// EUDAQ run end timestamp
  std::size_t eudaqRunEndTs;
  /// EPICS odd/even parity
  std::size_t epicsParity;
  /// EPICS PID
  std::size_t epicsPulseId;
  /// FACET-II DAQ number
  std::size_t epicsDAQNumber;
};

}  // namespace E320
