#include "TrackingPipeline/TrackFinding/SourceLinkGridConstructor.hpp"

#include "Acts/Surfaces/RectangleBounds.hpp"

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

SourceLinkGridConstructor::SourceLinkGridConstructor(const Config& config)
    : m_cfg(config) {};

std::unordered_map<Acts::GeometryIdentifier,
                   SourceLinkGridConstructor::GridType>
SourceLinkGridConstructor::constructGrid(
    const Acts::GeometryContext& gctx,
    std::vector<Acts::SourceLink> sourceLinks) {
  std::unordered_map<Acts::GeometryIdentifier, GridType> lookupTable;

  // Contruct lookup table
  for (const auto& [id, surface] : m_cfg.layers) {
    // Get bounds to construct the accumulator grid
    auto bounds =
        dynamic_cast<const Acts::RectangleBounds*>(&surface->bounds());

    if (bounds == nullptr) {
      throw std::invalid_argument("Only rectangle bounds supported");
    }
    if (surface->type() != Acts::Surface::SurfaceType::Plane) {
      throw std::invalid_argument("Only plane surfaces supported");
    }

    // Initialize the accumulator grid
    auto halfX = bounds->halfLengthX();
    auto halfY = bounds->halfLengthY();

    AxisType xAxis(-halfX, halfX, m_cfg.bins.first);
    AxisType yAxis(-halfY, halfY, m_cfg.bins.second);

    GridType grid(std::make_tuple(xAxis, yAxis));
    lookupTable.insert({id, grid});
  }
  // Fill the grid with source links
  for (auto& sl : sourceLinks) {
    auto ssl = sl.get<SimpleSourceLink>();
    auto sourceLinkId = ssl.geometryId().sensitive();

    Acts::GeometryIdentifier layerId;
    int layer = static_cast<int>(sourceLinkId / 10 - 1);
    layerId.setSensitive(layer + 1);

    // Convert hit from Detector surface to global
    Acts::Vector3 globalPos = m_cfg.surfaceAccessor(sl)->localToGlobal(
        gctx, ssl.parameters(), Acts::Vector3{0, 1, 0});

    // Convert hit from global to Layer surface
    Acts::Vector2 localPos =
        m_cfg.layers.at(layerId)
            ->globalToLocal(gctx, globalPos, Acts::Vector3(0, 1, 0), 1e-1)
            .value();

    auto bin = lookupTable.at(layerId).localBinsFromPosition(localPos);

    lookupTable.at(layerId).atLocalBins(bin).push_back(sl);
  }

  return lookupTable;
};
