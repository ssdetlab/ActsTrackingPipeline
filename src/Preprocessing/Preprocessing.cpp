#include "TrackingPipeline/Preprocessing/Preprocessing.hpp"

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "TrackingPipeline/Preprocessing/DetectorEvent.hpp"
#include "TrackingPipeline/Preprocessing/AnalysisFunctions.hpp"
#include "TrackingPipeline/Preprocessing/PairHash.hpp"
#include "TrackingPipeline/from_eudaq/DetectorEvent.h"
// ROOT
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"

namespace fs = std::filesystem;

namespace TrackingPipeline::Preprocessing {

void runPreprocessing(const PreprocessingConfig& cfg) {
  // 1. Build input TChain from all ROOT files in the specified directories
  TChain inChain(cfg.inputTreeName.c_str(), cfg.inputTreeName.c_str());

  for (const auto& inPath : cfg.inputDirs) {
    for (const auto& entry : fs::directory_iterator(inPath)) {
      if (!entry.is_regular_file() || entry.path().extension() != ".root") {
        continue;
      }
      inChain.Add(entry.path().c_str());
    }
  }

  const auto nEntries = static_cast<std::size_t>(inChain.GetEntries());
  if (nEntries <= cfg.skipEntries) {
    std::cerr << "Preprocessing: no entries to process (entries="
              << nEntries << ", skip=" << cfg.skipEntries << ")\n";
    return;
  }

  const double normEntries = static_cast<double>(nEntries - cfg.skipEntries);

  // 2. Set up input branch (original detector_event type)
  detector_event *detEvent = nullptr;
  inChain.SetBranchAddress(cfg.inputBranchName.c_str(), &detEvent);

  // 3. Set up output file and tree with ApollonIo::DetectorEvent
  TFile outFile(cfg.outputFile.c_str(), "RECREATE");
  TTree outTree(cfg.outputTreeName.c_str(), cfg.outputTreeName.c_str());

  ApollonIo::DetectorEvent outDetEvent;
  ULong64_t eventId = 0;
  outTree.Branch("event", &outDetEvent);
  outTree.Branch("eventId", &eventId);

  // 4. Histograms for raw and filtered occupancy (optional diagnostics)
  TH1I rawOcc("rawOcc", "rawOcc", 1024 * 512, 0, 1024 * 512);
  TH1I filteredOcc("filteredOcc", "filteredOcc", 1024 * 512, 0, 1024 * 512);

  // 5. First pass: accumulate per-chip hit statistics and raw occupancy
  std::unordered_map<int, double> norms;
  std::unordered_map<int, double> meanHitCounts;
  std::unordered_map<int, double> mean2HitCounts;
  std::unordered_map<int, std::unordered_map<Hit, double, PairHash>> totHitCounts;

  for (std::size_t i = cfg.skipEntries; i < nEntries; i++) {
    inChain.GetEntry(static_cast<Long64_t>(i));

    for (const auto &staveEv : detEvent->st_ev_buffer) {
      int staveId = 0;
      if (detEvent->run_number <= 169) {
        staveId = staveEv.stave_id;
      } else {
        staveId = (staveEv.stave_id == 0) ? 1 : 0;
      }

      for (const auto& chipEv : staveEv.ch_ev_buffer) {
        const int id = 100 * detEvent->run_number + 10 * staveId + chipEv.chip_id;

        norms[id] += 1.0;
        const double nhits = static_cast<double>(chipEv.hits.size());
        meanHitCounts[id] += nhits;
        mean2HitCounts[id] += nhits * nhits;

        for (const auto &[hitX, hitY] : chipEv.hits) {
          totHitCounts[id][{hitX, hitY}] += 1.0;
          rawOcc.Fill(512 * hitX + hitY);
        }
      }
    }
  }

  // 6. Identify noisy pixels per chip
  std::unordered_map<int, std::unordered_set<Hit, PairHash>> noisyPixels;

  for (const auto& [id, hitCounts] : totHitCounts) {
    const double mean = meanHitCounts[id] / norms[id];
    const double mean2 = mean2HitCounts[id] / norms[id];
    const double var = mean2 - mean * mean;  // variance (std^2)

    for (const auto &[hit, count] : hitCounts) {
      // Simple thresholds (could be made configurable)
      if (count > 1e4) {
        noisyPixels[id].insert(hit);
      }
      if (var > 0.0 && (count - mean) / var > 3.0) {
        noisyPixels[id].insert(hit);
      }
    }
  }

  // 7. Second pass: build filtered ApollonIo::DetectorEvent and filtered occupancy
  for (std::size_t i = cfg.skipEntries; i < nEntries; i++) {
    inChain.GetEntry(static_cast<Long64_t>(i));

    outDetEvent.trg_n = detEvent->trg_n;
    outDetEvent.ts_begin = detEvent->ts_begin;
    outDetEvent.ts_end = detEvent->ts_end;
    outDetEvent.run_number = detEvent->run_number;

    std::vector<ApollonIo::StaveEvent> outStaveEvents;
    outStaveEvents.reserve(detEvent->st_ev_buffer.size());

    // NOTE: this still uses the original single-stave logic
    const int requiredStaveId = 1;
    const int requiredChipSize = 5;

    for (const auto &staveEv : detEvent->st_ev_buffer) {
      int staveId = 0;
      if (detEvent->run_number <= 169) {
        staveId = staveEv.stave_id;
      } else {
        staveId = (staveEv.stave_id == 0) ? 1 : 0;
      }

      if (staveId != requiredStaveId) {
        continue;
      }
      if (static_cast<int>(staveEv.ch_ev_buffer.size()) != requiredChipSize) {
        continue;
      }

      ApollonIo::StaveEvent outStaveEvent;
      outStaveEvent.stave_id = staveId;

      std::vector<ApollonIo::ChipEvent> outChipEvents;
      outChipEvents.reserve(staveEv.ch_ev_buffer.size());

      for (const auto &chipEv : staveEv.ch_ev_buffer) {
        const int id = 100 * detEvent->run_number + 10 * staveId + chipEv.chip_id;
        ApollonIo::ChipEvent outChipEvent;

        // Copy flags and IDs
        outChipEvent.is_busy_violation = chipEv.is_busy_violation;
        outChipEvent.is_flushed_incomplete = chipEv.is_flushed_incomplete;
        outChipEvent.is_strobe_extended = chipEv.is_strobe_extended;
        outChipEvent.is_busy_transition = chipEv.is_busy_transition;

        outChipEvent.end_of_run = chipEv.end_of_run;
        outChipEvent.overflow = chipEv.overflow;
        outChipEvent.timeout = chipEv.timeout;
        outChipEvent.header_error = chipEv.header_error;
        outChipEvent.decoder_10b8b_error = chipEv.decoder_10b8b_error;
        outChipEvent.event_oversize_error = chipEv.event_oversize_error;

        outChipEvent.chip_id = chipEv.chip_id;
        outChipEvent.channel = chipEv.channel;

        // Filter noisy pixels and build hit set
        std::unordered_set<Hit, PairHash> hits;
        for (const auto &hit : chipEv.hits) {
          if (noisyPixels[id].contains(hit)) {
            continue;
          }
          hits.insert(hit);
          filteredOcc.Fill(512 * hit.first + hit.second);
        }
        if (hits.empty()) {
          continue;
        }

        // Cluster and store
        outChipEvent.hits = getClusters(hits);
        outChipEvents.push_back(std::move(outChipEvent));
      }

      if (static_cast<int>(outChipEvents.size()) != requiredChipSize) {
        continue;
      }

      outStaveEvent.ch_ev_buffer = std::move(outChipEvents);
      outStaveEvents.push_back(std::move(outStaveEvent));
    }

    if (outStaveEvents.empty()) {
      continue;
    }

    outDetEvent.st_ev_buffer = std::move(outStaveEvents);
    outTree.Fill();

    eventId++;
    if (i % 100 == 0) {
      std::cout << "Preprocessing: " << i << "/" << nEntries << "\n";
    }
  }

  // 8. Save diagnostic histograms
  const std::string baseDir =
    std::filesystem::path(cfg.outputFile).parent_path().string();

  TCanvas c1("c1", "c1", 800, 800);
  rawOcc.Draw("hist");
  c1.SaveAs((baseDir + "/raw.pdf").c_str());

  TCanvas c2("c2", "c2", 800, 800);
  filteredOcc.Draw("hist");
  c2.SaveAs((baseDir + "/filtered.pdf").c_str());

  // 9. Write tree and close file
  outTree.Write();
  outFile.Close();
}

}  // namespace TrackingPipeline::Preprocessing
