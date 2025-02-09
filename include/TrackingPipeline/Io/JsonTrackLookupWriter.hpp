#pragma once

#include <nlohmann/json.hpp>

#include "TrackingPipeline/Io/ITrackParamsLookupWriter.hpp"

/// @brief Json writer for track parameter lookup tables
///
/// This writer is used to write track parameter lookup tables
/// to a json file to be later used in track parameter estimation
/// for seeding
class JsonTrackLookupWriter final : public ITrackParamsLookupWriter {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Output file name
    std::string path;
  };

  /// Constructor
  ///
  /// @param config The configuration struct of the writer
  explicit JsonTrackLookupWriter(const Config& config);

  /// Write out track parameters lookup table
  ///
  /// @param lookup The lookup to write
  void writeLookup(const TrackLookup& lookup) const override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// The config of the writer
  Config m_cfg;
};
