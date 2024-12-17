#include "TrackingPipeline/Io/JsonTrackLookupReader.hpp"

#include "Acts/Plugins/Json/GridJsonConverter.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"

#include <fstream>

JsonTrackLookupReader::JsonTrackLookupReader(
    const Config& config) : m_cfg(config) {};

TrackLookup 
JsonTrackLookupReader::readLookup(const std::string& path) const {
    // Read the json file
    std::ifstream ifj(path);
    nlohmann::json jLookup;
    ifj >> jLookup;

    TrackLookup lookup;
    // Iterate over the json and deserialize the grids
    for (const auto& jGrid : jLookup) {
        Acts::GeometryIdentifier id(jGrid["geo_id"]);
    
        if (!m_cfg.refLayers.contains(id)) {
            throw std::invalid_argument("Geometry identifier not found");
        }
    
        const auto* refSurface = m_cfg.refLayers.at(id);
    
        // Get bounds to construct the lookup grid
        auto bounds =
            dynamic_cast<const Acts::RectangleBounds*>(&refSurface->bounds());
    
        if (bounds == nullptr) {
            throw std::invalid_argument("Only rectangle bounds supported");
        }
    
        // Axis is not deserilizable, so we need to recreate it
        auto halfX = bounds->halfLengthX();
        auto halfY = bounds->halfLengthY();
    
        TrackLookupAxisGen axisGen{
            {-halfX, halfX}, m_cfg.bins.first,
            {-halfY, halfY}, m_cfg.bins.second};
    
        // Deserialize the grid
        TrackLookupGrid grid =
            Acts::GridJsonConverter::fromJson<
                TrackLookupAxisGen, TrackLookupPair>(
                    jGrid["grid"], axisGen);

        lookup.try_emplace(id, std::move(grid));
    }

    return lookup;
};
