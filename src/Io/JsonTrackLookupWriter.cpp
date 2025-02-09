#include "TrackingPipeline/Io/JsonTrackLookupWriter.hpp"

#include "Acts/Plugins/Json/GridJsonConverter.hpp"

#include <fstream>

#include <nlohmann/json.hpp>

JsonTrackLookupWriter::JsonTrackLookupWriter(const Config& config)
    : m_cfg(config) {};

void JsonTrackLookupWriter::writeLookup(const TrackLookup& lookup) const {
  nlohmann::json jLookup;

  // Iterate over the lookup and serialize the grids
  for (const auto& [id, grid] : lookup) {
    nlohmann::json jGrid;
    jGrid["geo_id"] = id.value();
    jGrid["grid"] = Acts::GridJsonConverter::toJson(grid);

    jLookup.push_back(jGrid);
  }

  // Write the json file
  std::ofstream ofj(m_cfg.path, std::ios::out);
  ofj << std::setw(4) << jLookup << std::endl;
};
