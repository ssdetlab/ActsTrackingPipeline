#include "TrackingPipeline/Io/E320MagneticFieldWriter.hpp"

#include "Acts/Definitions/Units.hpp"

#include <memory>

#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/MagneticField/ConstantMagField.hpp"
#include "TrackingPipeline/MagneticField/IdealQuadrupoleMagField.hpp"

namespace E320 {

E320MagneticFieldWriter::E320MagneticFieldWriter(const Config& config,
                                                 Acts::Logging::Level level)
    : m_cfg(config) {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_file = new TFile(m_cfg.filePath.c_str(), "RECREATE");
  m_tree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());

  //------------------------------------------------------------------
  // Tree branches
  int bufSize = 32000;
  int splitLvl = 0;

  // Quad gradients T/m
  m_tree->Branch("quad1Grad", &m_quad1Grad, bufSize, splitLvl);
  m_tree->Branch("quad2Grad", &m_quad2Grad, bufSize, splitLvl);
  m_tree->Branch("quad3Grad", &m_quad3Grad, bufSize, splitLvl);

  /// Quad poitions mm
  m_tree->Branch("quad1CenterPrimary", &m_quad1CenterPrimary, bufSize,
                 splitLvl);
  m_tree->Branch("quad2CenterPrimary", &m_quad2CenterPrimary, bufSize,
                 splitLvl);
  m_tree->Branch("quad3CenterPrimary", &m_quad3CenterPrimary, bufSize,
                 splitLvl);

  m_tree->Branch("quad1CenterLong", &m_quad1CenterLong, bufSize, splitLvl);
  m_tree->Branch("quad2CenterLong", &m_quad2CenterLong, bufSize, splitLvl);
  m_tree->Branch("quad3CenterLong", &m_quad3CenterLong, bufSize, splitLvl);

  m_tree->Branch("quad1CenterShort", &m_quad1CenterShort, bufSize, splitLvl);
  m_tree->Branch("quad2CenterShort", &m_quad2CenterShort, bufSize, splitLvl);
  m_tree->Branch("quad3CenterShort", &m_quad3CenterShort, bufSize, splitLvl);

  // XCOR strength T
  m_tree->Branch("xCorrectorStrength", &m_xCorrectorStrength, bufSize,
                 splitLvl);

  // Dipole strength T
  m_tree->Branch("dipoleStrength", &m_dipoleStrength, bufSize, splitLvl);

  // Event number handle
  m_tree->Branch("eventId", &m_eventId, bufSize, splitLvl);
}

ProcessCode E320MagneticFieldWriter::finalize() {
  if (m_file != nullptr) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode E320MagneticFieldWriter::write(const AlgorithmContext& ctx) {
  std::lock_guard<std::mutex> lock(m_mutex);

  m_eventId = ctx.eventNumber;

  const auto& mctx = ctx.magFieldContext;
  const auto& goInst = *E320::GeometryOptions::instance();

  // Quad gradients T/m
  m_quad1Grad =
      goInst.quad1Gradient * Acts::UnitConstants::m / Acts::UnitConstants::T;
  m_quad2Grad =
      goInst.quad2Gradient * Acts::UnitConstants::m / Acts::UnitConstants::T;
  m_quad3Grad =
      goInst.quad3Gradient * Acts::UnitConstants::m / Acts::UnitConstants::T;

  // Quad poitions mm
  m_quad1CenterPrimary = goInst.quad1CenterPrimary / Acts::UnitConstants::mm;
  m_quad2CenterPrimary = goInst.quad2CenterPrimary / Acts::UnitConstants::mm;
  m_quad3CenterPrimary = goInst.quad3CenterPrimary / Acts::UnitConstants::mm;

  m_quad1CenterLong = goInst.quad1CenterLong / Acts::UnitConstants::mm;
  m_quad2CenterLong = goInst.quad2CenterLong / Acts::UnitConstants::mm;
  m_quad3CenterLong = goInst.quad3CenterLong / Acts::UnitConstants::mm;

  m_quad1CenterShort = goInst.quad1CenterShort / Acts::UnitConstants::mm;
  m_quad2CenterShort = goInst.quad2CenterShort / Acts::UnitConstants::mm;
  m_quad3CenterShort = goInst.quad3CenterShort / Acts::UnitConstants::mm;

  // XCOR strength T
  m_xCorrectorStrength =
      goInst.xCorrectorFieldStrength / Acts::UnitConstants::T;

  // Dipole strength T
  m_dipoleStrength = goInst.dipoleFieldStrength / Acts::UnitConstants::T;

  if (mctx.hasValue()) {
    const auto& store = mctx.get<std::shared_ptr<MagneticFieldStore>&>();

    if (store->store.contains(goInst.quad1Id)) {
      m_quad1Grad = store->store.at(goInst.quad1Id)
                        .as<IdealQuadrupoleMagField::Cache>()
                        .m_gradient *
                    Acts::UnitConstants::m / Acts::UnitConstants::T;
    }
    if (store->store.contains(goInst.quad2Id)) {
      m_quad2Grad = store->store.at(goInst.quad2Id)
                        .as<IdealQuadrupoleMagField::Cache>()
                        .m_gradient *
                    Acts::UnitConstants::m / Acts::UnitConstants::T;
    }
    if (store->store.contains(goInst.quad3Id)) {
      m_quad3Grad = store->store.at(goInst.quad3Id)
                        .as<IdealQuadrupoleMagField::Cache>()
                        .m_gradient *
                    Acts::UnitConstants::m / Acts::UnitConstants::T;
    }
    if (store->store.contains(goInst.xCorrectorId)) {
      m_xCorrectorStrength = store->store.at(goInst.xCorrectorId)
                                 .as<ConstantMagField::Cache>()
                                 .m_field(goInst.xCorrectorDirIdx) /
                             Acts::UnitConstants::T;
    }
    if (store->store.contains(goInst.dipoleId)) {
      m_dipoleStrength = store->store.at(goInst.dipoleId)
                             .as<ConstantMagField::Cache>()
                             .m_field(goInst.dipoleDirIdx) /
                         Acts::UnitConstants::T;
    }
  }

  m_tree->Fill();

  // Return success flag
  return ProcessCode::SUCCESS;
}

}  // namespace E320
