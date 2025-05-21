#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"

class E320DipoleTrackLookupProvider {
 public:
  struct Config {
    /// Reference surface
    const Acts::Surface* referenceSurface;

    double layerPosition;

    double dipoleSize;
    double dipolePosition;
    double dipoleAmplidute;

    double correctorSize;
    double correctorPosition;
    double correctorAmplidute;
  };

  E320DipoleTrackLookupProvider(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);

  /// Lookup the track parameters at a given position
  std::pair<Acts::CurvilinearTrackParameters, Acts::CurvilinearTrackParameters>
  lookup(const Acts::GeometryContext& gctx,
         const Acts::SourceLink& pivot) const;

  /// @brief Readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  const Acts::Logger& logger() const { return *m_logger; }

 private:
  double m_layerDipoleDistance;
  double m_layerCorrectorDistance;
  Acts::BoundSquareMatrix m_cov;

  std::unique_ptr<const Acts::Logger> m_logger;

  /// Configuration
  Config m_cfg;
};
