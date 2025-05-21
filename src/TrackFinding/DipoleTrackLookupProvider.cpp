#include "TrackingPipeline/TrackFinding/DipoleTrackLookupProvider.hpp"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <stdexcept>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

using namespace Acts::UnitLiterals;

E320DipoleTrackLookupProvider::E320DipoleTrackLookupProvider(
    const Config& cfg, Acts::Logging::Level level)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("DipoleTrackLookupProvider", level)) {
  if (m_cfg.referenceSurface == nullptr) {
    throw std::runtime_error("Reference surface is not initialized");
  }
  Acts::BoundVector ipStdDev;
  ipStdDev[Acts::eBoundLoc0] = 100_um;
  ipStdDev[Acts::eBoundLoc1] = 100_um;
  ipStdDev[Acts::eBoundTime] = 25_ns;
  ipStdDev[Acts::eBoundPhi] = 2_degree;
  ipStdDev[Acts::eBoundTheta] = 2_degree;
  ipStdDev[Acts::eBoundQOverP] = 1 / 100_GeV;
  m_cov = ipStdDev.cwiseProduct(ipStdDev).asDiagonal();

  m_layerDipoleDistance = (m_cfg.layerPosition - m_cfg.dipolePosition);
  m_layerCorrectorDistance = (m_cfg.layerPosition - m_cfg.correctorPosition);
  ACTS_INFO("initialized with reference layer "
            << m_cfg.referenceSurface->geometryId());
  ACTS_INFO("Layer-dipole distance " << m_layerDipoleDistance << " mm");
  ACTS_INFO("Layer-corrector distance " << m_layerCorrectorDistance << " mm");
  ACTS_INFO("Dipole amplitude " << m_cfg.dipoleAmplidute << " T");
  ACTS_INFO("Dipole size " << m_cfg.dipoleSize << " m");
}

std::pair<Acts::CurvilinearTrackParameters, Acts::CurvilinearTrackParameters>
E320DipoleTrackLookupProvider::lookup(const Acts::GeometryContext& gctx,
                                      const Acts::SourceLink& pivot) const {
  const auto& ssl = pivot.get<SimpleSourceLink>();

  ACTS_VERBOSE("Analysing measurement at GeoID " << ssl.geometryId());
  if (ssl.geometryId() != m_cfg.referenceSurface->geometryId()) {
    throw std::runtime_error("Pivot is not on the reference surface");
  }

  Acts::Vector3 globalPos = m_cfg.referenceSurface->localToGlobal(
      gctx, ssl.parameters(), Acts::Vector3::UnitY());

  double phi = std::atan(globalPos.x() / m_layerCorrectorDistance);
  double theta = std::atan(globalPos.z() / m_layerDipoleDistance);
  double P = std::abs(0.3 * m_cfg.dipoleAmplidute * m_cfg.dipoleSize /
                      (std::sin(theta) * std::cos(phi)));

  ACTS_VERBOSE("Global measurement position " << globalPos.transpose());
  ACTS_VERBOSE("Reference layer momentum theta " << theta);
  ACTS_VERBOSE("Momentum magnitude " << P);

  /*Acts::CurvilinearTrackParameters ipPars(*/
  /*    Acts::Vector4(0, 0, 0, 0), Acts::Vector3(0, 1, 0), 1_e / P, m_cov,*/
  /*    Acts::ParticleHypothesis::electron());*/
  /*Acts::CurvilinearTrackParameters ipPars(*/
  /*    Acts::Vector4(globalPos.x() - 126_mm * std::cos(theta) *
   * std::sin(phi),*/
  /*                  globalPos.y() - 126_mm * std::cos(theta) *
   * std::cos(phi),*/
  /*                  globalPos.z() - 126_mm * std::sin(theta), 0),*/
  /*    Acts::Vector3(std::cos(theta) * std::sin(phi),*/
  /*                  std::cos(theta) * std::cos(phi), std::sin(theta)),*/
  /*    1_e / P, m_cov, Acts::ParticleHypothesis::electron());*/
  Acts::CurvilinearTrackParameters ipPars(
      Acts::Vector4(globalPos.x() - 1_mm * std::cos(theta) * std::sin(phi),
                    globalPos.y() - 1_mm * std::cos(theta) * std::cos(phi),
                    globalPos.z() - 1_mm * std::sin(theta), 0),
      Acts::Vector3(std::cos(theta) * std::sin(phi),
                    std::cos(theta) * std::cos(phi), std::sin(theta)),
      1_e / P, m_cov, Acts::ParticleHypothesis::electron());
  Acts::CurvilinearTrackParameters refPars(
      Acts::Vector4(globalPos.x(), globalPos.y(), globalPos.z(), 0),
      Acts::Vector3(std::cos(theta) * std::sin(phi),
                    std::cos(theta) * std::cos(phi), std::sin(theta)),
      1_e / P, m_cov, Acts::ParticleHypothesis::electron());

  return {ipPars, refPars};
}
