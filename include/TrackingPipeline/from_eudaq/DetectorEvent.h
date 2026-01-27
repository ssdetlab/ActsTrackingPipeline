#ifndef DetectorEvent_h
#define DetectorEvent_h

#include <TObject.h>
#include <cstdint>
#include <vector>

#include "ChipRegs.h"
#include "MOSAICRegs.h"

class chip_event {
public:
  /// Hit map of the chip in the event
  std::vector<std::pair<std::uint16_t, std::uint16_t>> hits;

  /// ALPIDE readout flags
  bool is_busy_violation;
  bool is_flushed_incomplete;
  bool is_strobe_extended;
  bool is_busy_transition;

  /// MOSAIC readout flags
  bool end_of_run;
  bool overflow;
  bool timeout;
  bool header_error;
  bool decoder_10b8b_error;
  bool event_oversize_error;

  /// Chip identificators
  std::uint8_t chip_id;
  std::uint8_t channel;

  /// Temperature
  float temperature;
  float delta_temperature;

  ClassDef(chip_event, 5);
};

class stave_event {
public:
  /// Chip event storage
  std::vector<chip_event> ch_ev_buffer;

  /// Event MOSAIC trigger ID
  std::uint64_t trg_n;

  /// Stave identificator
  std::uint8_t stave_id;

  /// UNIX timestamp
  std::uint64_t ts_begin;
  std::uint64_t ts_end;

  ClassDef(stave_event, 4);
};

class detector_event {
public:
  /// Stave event storage
  std::vector<stave_event> st_ev_buffer;

  /// Event MOSAIC trigger ID
  std::uint64_t trg_n;

  /// UNIX timestamp
  std::uint64_t ts_begin;
  std::uint64_t ts_end;

  /// Run number
  std::uint32_t run_number;

  ClassDef(detector_event, 7);
};

class chip_run_meta {
public:
  /// Chip identificator
  std::uint8_t chip_id;

  /// ALPIDE conf
  ALPIDERegs chip_regs;

  /// Number of received triggers
  std::uint16_t n_trigs;

  /// Number of generated strobes
  std::uint16_t n_strobes;

  /// Number of matrix readouts
  std::uint16_t n_matrix_readouts;

  /// Number of sent frames
  std::uint16_t n_frames;

  ClassDef(chip_run_meta, 2);
};

class stave_run_meta {
public:
  /// Stave identificator
  std::uint8_t stave_id;

  /// Frames with data
  std::uint32_t non_empty_frames;

  /// Number of triggers received
  /// By MOSAIC
  std::uint32_t n_trigs_mosaic;

  /// MOSAIC conf
  MOSAICRegs mosaic_regs;

  /// ALPIDE meta data
  std::vector<chip_run_meta> chips_rmd;

  ClassDef(stave_run_meta, 2);
};

class run_meta {
public:
  using time_server_md =
      std::tuple<std::uint64_t, std::uint64_t, std::uint64_t>;

  /// Current run number
  std::uint32_t run_number;

  /// UNIX timestamp in ms
  std::uint64_t run_start;
  std::uint64_t run_end;

  /// Stave meta data
  std::vector<stave_run_meta> st_rmd;

  /// Time server meta data
  std::vector<time_server_md> ts_md;

  ClassDef(run_meta, 2);
};

#endif