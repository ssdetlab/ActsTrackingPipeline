#include "TrackingPipeline/TrackFinding/TrackLookupProvider.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/EventData/ParticleHypothesis.hpp>

#include <array>
#include <optional>
#include <stdexcept>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

TrackLookupProvider::TrackLookupProvider(const Config& config)
    : m_cfg(std::move(config)),
      m_lookup(std::make_shared<TrackLookup>(
          m_cfg.trackLookupReader->readLookup(m_cfg.lookupPath))) {}

// TODO: Prioritize direction along x
std::optional<TrackLookupGrid::index_t> TrackLookupProvider::findClosestFilled(
    const TrackLookupGrid& grid, TrackLookupGrid::index_t bin,
    const Acts::Vector2& dir) const {
  std::vector<TrackLookupGrid::index_t> neighbours;

  // Origin
  if (dir == Acts::Vector2{0, 0}) {
    neighbours = {
        {bin.at(0) - 1, bin.at(1) + 1}, {bin.at(0), bin.at(1) + 1},
        {bin.at(0) + 1, bin.at(1) + 1}, {bin.at(0) - 1, bin.at(1)},
        {bin.at(0) + 1, bin.at(1)},     {bin.at(0) - 1, bin.at(1) - 1},
        {bin.at(0), bin.at(1) - 1},     {bin.at(0) + 1, bin.at(1) - 1}};
  } else if (dir == Acts::Vector2{-1, 1}) {
    neighbours = {{bin.at(0) - 1, bin.at(1) + 1},
                  {bin.at(0), bin.at(1) + 1},
                  {bin.at(0) - 1, bin.at(1)}};
  }
  if (dir == Acts::Vector2{0, 1}) {
    neighbours = {{bin.at(0), bin.at(1) + 1}};
  }
  if (dir == Acts::Vector2{1, 1}) {
    neighbours = {{bin.at(0), bin.at(1) + 1},
                  {bin.at(0) + 1, bin.at(1) + 1},
                  {bin.at(0) + 1, bin.at(1)}};
  }
  if (dir == Acts::Vector2{1, 0}) {
    neighbours = {{bin.at(0) + 1, bin.at(1)}};
  }
  if (dir == Acts::Vector2{1, -1}) {
    neighbours = {{bin.at(0) + 1, bin.at(1)},
                  {bin.at(0), bin.at(1) - 1},
                  {bin.at(0) + 1, bin.at(1) - 1}};
  }
  if (dir == Acts::Vector2{0, -1}) {
    neighbours = {{bin.at(0), bin.at(1) - 1}};
  }
  if (dir == Acts::Vector2{-1, -1}) {
    neighbours = {{bin.at(0) - 1, bin.at(1)},
                  {bin.at(0) - 1, bin.at(1) - 1},
                  {bin.at(0), bin.at(1) - 1}};
  }

  for (const auto& n : neighbours) {
    if (grid.atLocalBins(n).first != nullptr) {
      return std::optional<TrackLookupGrid::index_t>(n);
    }
  }
  for (const auto& n : neighbours) {
    Acts::Vector2 newDir{static_cast<std::int32_t>(n.at(0)) -
                             static_cast<std::int32_t>(bin.at(0)),
                         static_cast<std::int32_t>(n.at(1)) -
                             static_cast<std::int32_t>(bin.at(1))};
    auto res = findClosestFilled(grid, n, newDir);
    if (res.has_value()) {
      return res;
    }
  }
  return std::nullopt;
}

std::pair<Acts::CurvilinearTrackParameters, Acts::CurvilinearTrackParameters>
TrackLookupProvider::lookup(const Acts::GeometryContext& gctx,
                            const Acts::SourceLink& pivot) const {
  auto ssl = pivot.get<SimpleSourceLink>();

  Acts::Vector2 localPos = ssl.parametersLoc();
  Acts::GeometryIdentifier geoId = ssl.geometryId();

  auto bin = m_lookup->at(geoId).localBinsFromPosition(localPos);

  if (m_lookup->at(geoId).atLocalBins(bin).first == nullptr) {
    try {
      auto res = findClosestFilled(m_lookup->at(geoId), bin, {0, 0});
      bin = res.value();
    } catch (const std::out_of_range& err) {
      throw std::runtime_error("Cannot find filled bin in a grid");
    }
  }
  return {*m_lookup->at(geoId).atLocalBins(bin).first,
          *m_lookup->at(geoId).atLocalBins(bin).second};
}
