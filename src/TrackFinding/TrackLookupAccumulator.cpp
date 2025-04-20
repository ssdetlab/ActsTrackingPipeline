#include "TrackingPipeline/TrackFinding/TrackLookupAccumulator.hpp"

TrackLookupAccumulator::TrackLookupAccumulator(TrackLookupGrid grid)
    : m_grid(std::move(grid)) {}

void TrackLookupAccumulator::addTrack(
    const Acts::CurvilinearTrackParameters& ipTrackParameters,
    const Acts::CurvilinearTrackParameters& refTrackParameters,
    const Acts::Vector2& position) {
  std::lock_guard<std::mutex> lock(m_gridMutex);

  auto bin = m_grid.localBinsFromPosition(position);

  if (m_countGrid[bin] == 0) {
    m_grid.atLocalBins(bin).first =
        std::make_shared<Acts::CurvilinearTrackParameters>(ipTrackParameters);
    m_grid.atLocalBins(bin).second =
        std::make_shared<Acts::CurvilinearTrackParameters>(refTrackParameters);

    m_countGrid.at(bin)++;
    return;
  }

  *m_grid.atLocalBins(bin).first =
      addTrackParameters(*m_grid.atLocalBins(bin).first, ipTrackParameters);
  *m_grid.atLocalBins(bin).second =
      addTrackParameters(*m_grid.atLocalBins(bin).second, refTrackParameters);
  m_countGrid.at(bin)++;
}

TrackLookupGrid TrackLookupAccumulator::finalizeLookup() {
  auto meanTrack = [&](const Acts::CurvilinearTrackParameters& track,
                       std::size_t count) {
    return Acts::CurvilinearTrackParameters(
        track.fourPosition() / count, track.momentum().normalized(),
        count * track.charge() / track.momentum().norm(), track.covariance(),
        track.particleHypothesis());
  };

  int totCount = 0;
  for (auto [bin, count] : m_countGrid) {
    if (count == 0) {
      continue;
    }
    *m_grid.atLocalBins(bin).first =
        meanTrack(*m_grid.atLocalBins(bin).first, count);
    *m_grid.atLocalBins(bin).second =
        meanTrack(*m_grid.atLocalBins(bin).second, count);
    totCount += count;
  }

  return m_grid;
}

Acts::CurvilinearTrackParameters TrackLookupAccumulator::addTrackParameters(
    const Acts::CurvilinearTrackParameters& a,
    const Acts::CurvilinearTrackParameters& b) {
  if (a.particleHypothesis() != b.particleHypothesis()) {
    throw std::invalid_argument(
        "Cannot accumulate track parameters with different particle "
        "hypotheses");
  }
  if (a.charge() != b.charge()) {
    throw std::invalid_argument(
        "Cannot accumulate track parameters with different charges");
  }

  Acts::Vector3 momentum = a.momentum() + b.momentum();
  Acts::Vector4 fourPosition = a.fourPosition() + b.fourPosition();
  return Acts::CurvilinearTrackParameters(
      fourPosition, momentum.normalized(), a.charge() / momentum.norm(),
      a.covariance(), a.particleHypothesis());
}
