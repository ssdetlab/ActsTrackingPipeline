#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"

class E320DipoleTrackLookupProvider {
 public:
  struct Config {
    /// Reference surface
    const Acts::Surface* referenceSurface;

    /// Reference surface z-position in mm
    double layerPosition;

    /// Detector tilt in yz-plane as a whole in rad
    double detectorYZTilt;

    /// Dipole z-length in m
    double dipoleSize;
    /// Dipole z-position in mm
    double dipolePosition;
    /// Dipole field amplitude in T
    double dipoleAmplidute;

    /// Corrector z-length in m
    double correctorSize;
    /// Corrector z-position in mm
    double correctorPosition;
    /// Corrector field amplitude in T
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
