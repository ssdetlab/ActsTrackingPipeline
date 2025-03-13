#pragma once

#include <nlohmann/json.hpp>

#include "TrackingPipeline/Io/ITrackParamsLookupReader.hpp"

/// @brief Json reader for track parameter lookup tables
///
/// This reader is used to read track parameter lookup tables
/// from a json file to be later used in track parameter estimation
/// for seeding
class JsonTrackLookupReader final : public ITrackParamsLookupReader {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Reference tracking layers
    std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
        refLayers;
    /// Binning of the grid to be emposed
    /// onto the reference layers
    std::pair<std::size_t, std::size_t> bins;
  };

  explicit JsonTrackLookupReader(const Config& config);

  /// @brief Read the lookup from a json file
  ///
  /// @param path path to the json file
  ///
  /// @return lookup table for track parameter estimation
  TrackLookup readLookup(const std::string& path) const override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration struct
  Config m_cfg;
};
