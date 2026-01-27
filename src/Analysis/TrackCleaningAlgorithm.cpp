#include "TrackingPipeline/Analysis/TrackCleaningAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmRegistry.hpp"

#include <toml.hpp>
#include <algorithm>
#include <unordered_map>

//---------------------------------------------------------------------
// Constructor
//---------------------------------------------------------------------

TrackCleaningAlgorithm::TrackCleaningAlgorithm(const Config& cfg,
                                               Acts::Logging::Level lvl)
  : IAlgorithm("TrackCleaningAlgorithm", lvl)
  , m_cfg(cfg) {

  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("TrackCleaningAlgorithm: outputTracks empty");
  }

  if (!m_cfg.inputCleaningTracks.empty()) {
    m_inputCleaning.initialize(m_cfg.inputCleaningTracks);
  }
  if (!m_cfg.inputActsTracks.empty()) {
    m_inputActs.initialize(m_cfg.inputActsTracks);
  }

  if (m_cfg.inputCleaningTracks.empty() && m_cfg.inputActsTracks.empty()) {
    throw std::invalid_argument(
        "TrackCleaningAlgorithm: at least one of inputCleaningTracks or "
        "inputActsTracks must be set");
  }

  m_outputTracks.initialize(m_cfg.outputTracks);
}

//---------------------------------------------------------------------
// Helper: convert Acts Tracks -> CleaningTracks
//---------------------------------------------------------------------

namespace {

// Build CleaningTracks (TrackDescriptor vector) from Acts Tracks,
// using the same chi2Smoothed and hit extraction logic as RootTrackWriter.
static CleaningTracks
makeCleaningTracksFromActs(const Tracks& inputTracks,
                           const AlgorithmContext& ctx) {
  CleaningTracks out;
  out.reserve(inputTracks.tracks.size());

  // Use eventNumber from context as eventId (you can refine this if needed)
  std::size_t eventId = ctx.eventNumber;

  // Loop over fitted tracks
  for (std::size_t tid = 0; tid < inputTracks.tracks.size(); ++tid) {
    const auto& track = inputTracks.tracks.getTrack(tid);

    TrackDescriptor desc;

    desc.eventId = eventId;
    desc.trackId = static_cast<std::size_t>(inputTracks.trackIds.at(tid));

    // You can fill these later if you have them
    desc.pdgId  = 0;
    desc.charge = 0;

    // ndf from Acts track
    desc.ndf = track.nDoF();

    // chi2Smoothed and global hit positions, following RootTrackWriter
    double chi2Smoothed = 0.0;
    std::vector<Acts::Vector3> hitsGlobal;
    hitsGlobal.reserve(track.nTrackStates());

    for (const auto& state : track.trackStatesReversed()) {
      // Skip states without measurements / projector
      if (!state.hasProjector()) {
        continue;
      }

      // Measurement in local coordinates
      Acts::Vector2 hit = state.effectiveCalibrated();

      // Global position of measurement hit
      Acts::Vector3 hitGlobal =
          state.referenceSurface().localToGlobal(
              ctx.geoContext, hit, Acts::Vector3(1., 0., 0.));
      hitsGlobal.push_back(hitGlobal);

      // We follow the same chi2Smoothed definition as RootTrackWriter
      if (state.hasSmoothed()) {
        Acts::Vector2 smoothedHit =
            state.effectiveProjector() * state.smoothed();

        Acts::Vector2 smoothedResidual = hit - smoothedHit;

        Acts::SquareMatrix2 measurementCov =
            state.effectiveCalibratedCovariance();

        Acts::SquareMatrix2 smoothedCov =
            state.effectiveProjector() * state.smoothedCovariance() *
                state.effectiveProjector().transpose() -
            measurementCov;

        Acts::Vector2 smoothedDiag =
            smoothedCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

        Acts::Vector2 smoothedPull =
            smoothedDiag.cwiseProduct(smoothedResidual);

        chi2Smoothed += smoothedPull.dot(smoothedPull);
      }
    }

    desc.chi2Smoothed    = chi2Smoothed;
    desc.trackHitsGlobal = std::move(hitsGlobal);

    out.push_back(std::move(desc));
  }

  return out;
}

} // namespace

//---------------------------------------------------------------------
// execute()
//---------------------------------------------------------------------

ProcessCode TrackCleaningAlgorithm::execute(const AlgorithmContext& ctx) const {
  CleaningContainer workingTracks;

  // Case 1: cleaning-only mode – CleaningTracks directly from reader
  if (!m_cfg.inputCleaningTracks.empty()) {
    const auto& in = m_inputCleaning(ctx);
    if (!in.empty()) {
      workingTracks = in;
    }
  }

  // Case 2: full run – Acts Tracks, convert on the fly
  if (workingTracks.empty() && !m_cfg.inputActsTracks.empty()) {
    const auto& actsIn = m_inputActs(ctx);
    workingTracks = makeCleaningTracksFromActs(actsIn, ctx);
  }

  if (workingTracks.empty()) {
    m_outputTracks(ctx, CleaningContainer{});
    return ProcessCode::SUCCESS;
  }

  // Group by eventId and sort by chi2Smoothed
  std::unordered_map<std::size_t,
      std::vector<const CleaningContainer::value_type*>> perEvent;
  perEvent.reserve(workingTracks.size());

  for (const auto& t : workingTracks) {
    perEvent[t.eventId].push_back(&t);
  }

  for (auto& [id, vec] : perEvent) {
    std::sort(vec.begin(), vec.end(),
              [](const auto* a, const auto* b) {
                return a->chi2Smoothed < b->chi2Smoothed;
              });
  }

  CleaningContainer cleaned;
  cleaned.reserve(workingTracks.size());

  // For each event: keep best chi2 track that doesn't share hits
  for (auto& [id, vec] : perEvent) {
    std::vector<const CleaningContainer::value_type*> unique;
    unique.reserve(vec.size());

    for (auto* ta : vec) {
      bool shared = false;
      for (auto* tb : unique) {
        const auto& uniqueHits = tb->trackHitsGlobal;
        for (const auto& hit : ta->trackHitsGlobal) {
          if (std::find(uniqueHits.begin(), uniqueHits.end(), hit) !=
              uniqueHits.end()) {
            shared = true;
            break;
          }
        }
        if (shared) break;
      }
      if (!shared) {
        unique.push_back(ta);
      }
    }

    for (auto* t : unique) {
      cleaned.push_back(*t);
    }
  }

  m_outputTracks(ctx, std::move(cleaned));
  return ProcessCode::SUCCESS;
}

//---------------------------------------------------------------------
// Registrar
//---------------------------------------------------------------------

namespace {

struct TrackCleaningAlgorithmRegistrar {
  TrackCleaningAlgorithmRegistrar() {
    using namespace TrackingPipeline;

    AlgorithmRegistry::instance().registerBuilder(
      "TrackCleaningAlgorithm",
      [](const toml::value& section,
         Acts::Logging::Level logLevel) -> AlgorithmPtr {

        TrackCleaningAlgorithm::Config cfg;
        cfg.inputCleaningTracks =
            toml::find_or<std::string>(section, "inputCleaningTracks", "");
        cfg.inputActsTracks =
            toml::find_or<std::string>(section, "inputActsTracks", "");
        cfg.outputTracks =
            toml::find<std::string>(section, "outputTracks");

        return std::make_shared<TrackCleaningAlgorithm>(cfg, logLevel);
      });
  }
} _TrackCleaningAlgorithmRegistrar;

}  // namespace
