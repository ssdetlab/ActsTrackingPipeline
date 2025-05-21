#ifndef DetectorEvent_h
#define DetectorEvent_h

#include <cstddef>
#include <set>
#include <tuple>

#include <TObject.h>

#include "AlpideRegisters.h"
#include "EpicsFrame.h"
#include "MosaicRegisters.h"

namespace E320Io {

class ChipEvent {
 public:
  struct Cluster {
    double xCenter;
    double yCenter;

    std::size_t dx;
    std::size_t dy;
    std::size_t cl_size;

    bool operator<(const Cluster& cl) const {
      return std::tie(xCenter, yCenter, dx, dy, cl_size) <
             std::tie(cl.xCenter, cl.yCenter, cl.dx, cl.dx, cl.cl_size);
    }
  };

  // hit map of the chip in the event
  /*std::set<std::tuple<double, double, std::size_t, std::size_t, std::size_t>>*/
  /*    hits;*/
  std::set<Cluster> hits;

  // ALPIDE readout flags
  bool is_busy_violation;
  bool is_flushed_incomplete;
  bool is_strobe_extended;
  bool is_busy_transition;

  // MOSAIC readout flags
  bool end_of_run;
  bool overflow;
  bool timeout;
  bool header_error;
  bool decoder_10b8b_error;
  bool event_oversize_error;

  // chip identificators
  std::uint8_t chip_id;
  std::uint8_t channel;

  ClassDef(ChipEvent, 3);
};

class StaveEvent {
 public:
  // chip event storage
  std::vector<ChipEvent> ch_ev_buffer;

  // stave identificator
  std::uint8_t stave_id;

  ClassDef(StaveEvent, 3);
};

class DetectorEvent {
 public:
  // stave event storage
  std::vector<StaveEvent> st_ev_buffer;

  // Epics frame
  EpicsFrame epics_frame;

  // event MOSAIC trigger ID
  std::uint64_t trg_n;

  // UNIX timestamp in ms
  std::uint64_t ts_begin;
  std::uint64_t ts_end;

  // Run number
  std::uint32_t run_number;

  ClassDef(DetectorEvent, 6);
};

class ChipRunMeta {
 public:
  // chip identificator
  std::uint8_t chip_id;

  // ALPIDE conf
  ALPIDERegs chip_regs;

  // number of received triggers
  std::uint16_t n_trigs;

  // number of generated strobes
  std::uint16_t n_strobes;

  // number of matrix readouts
  std::uint16_t n_matrix_readouts;

  // number of sent frames
  std::uint16_t n_frames;

  ClassDef(ChipRunMeta, 2);
};

class StaveRunMeta {
 public:
  // stave identificator
  std::uint8_t stave_id;

  // frames with data
  std::uint32_t non_empty_frames;

  // number of triggers received
  // by MOSAIC
  std::uint32_t n_trigs_mosaic;

  // MOSAIC conf
  MOSAICRegs mosaic_regs;

  // ALPIDE meta data
  std::vector<ChipRunMeta> chips_rmd;

  ClassDef(StaveRunMeta, 2);
};

class RunMeta {
 public:
  using time_server_md =
      std::tuple<std::uint64_t, std::uint64_t, std::uint64_t>;

  // current run number
  std::uint32_t run_number;

  // UNIX timestamp in ms
  std::uint64_t run_start;
  std::uint64_t run_end;

  // Stave meta data
  std::vector<StaveRunMeta> st_rmd;

  // time server meta data
  std::vector<time_server_md> ts_md;

  ClassDef(RunMeta, 2);
};

}  // namespace E320Io

#endif
