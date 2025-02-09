#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"

#include "TrackingPipeline/Io/ITrackParamsLookupReader.hpp"

/// @brief Class that provides estimates of IP and
/// reference layer track parameters given pivot source link
///
/// Class reads estimated track parameters with the provided reader
/// and provides estimates of the IP and reference layer track
/// parameters given the pivot source link
class TrackLookupProvider {
 public:
  struct Config {
    /// Lookup reader
    std::shared_ptr<ITrackParamsLookupReader> trackLookupReader;
    /// Lookup path
    std::string lookupPath;
  };

  TrackLookupProvider(const Config& config);

  /// Lookup the track parameters at a given position
  std::pair<Acts::CurvilinearTrackParameters, Acts::CurvilinearTrackParameters>
  lookup(const Acts::GeometryContext& gctx,
         const Acts::SourceLink& pivot) const;

  /// @brief Readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  std::optional<TrackLookupGrid::index_t> findClosestFilled(
      const TrackLookupGrid& grid, TrackLookupGrid::index_t bin,
      const Acts::Vector2& dir) const;

  /// Configuration
  Config m_cfg;

  std::shared_ptr<TrackLookup> m_lookup;
};
