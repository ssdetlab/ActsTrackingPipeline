#include "TrackingPipeline/Io/E320MagneticFieldParametersProvider.hpp"

#include "Acts/Definitions/Units.hpp"

#include <cstddef>
#include <memory>

#include "TrackingPipeline/Geometry/E320GeometryOptions.hpp"
#include "TrackingPipeline/MagneticField/ConstantMagField.hpp"
#include "TrackingPipeline/MagneticField/IdealQuadrupoleMagField.hpp"

namespace E320 {

E320MagneticFieldParametersProvider::E320MagneticFieldParametersProvider(
    const Config& config)
    : m_cfg(config) {
  if (m_cfg.filePaths.empty()) {
    throw std::invalid_argument("Missing input filenames");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  if (m_cfg.filePaths.size() == 1) {
    m_file = new TFile(m_cfg.filePaths.at(0).c_str());
    m_tree = m_file->Get<TTree>(m_cfg.treeName.c_str());
  } else {
    m_chainOwner = new TChain(m_cfg.treeName.c_str());
    // Add the files to the chain
    for (const auto& path : m_cfg.filePaths) {
      m_chainOwner->Add(path.c_str());
    }
    m_tree = dynamic_cast<TTree*>(m_chainOwner);
  }

  // Quad gradients T/m
  m_tree->SetBranchAddress("quad1Grad", &m_quad1Grad);
  m_tree->SetBranchAddress("quad2Grad", &m_quad2Grad);
  m_tree->SetBranchAddress("quad3Grad", &m_quad3Grad);

  // XCOR strength T
  m_tree->SetBranchAddress("xCorrectorStrength", &m_xCorrectorStrength);

  // Dipole strength T
  m_tree->SetBranchAddress("dipoleStrength", &m_dipoleStrength);

  // Event Id
  m_tree->SetBranchAddress("eventId", &m_eventId);

  // Geometry constraints instance
  const auto& goInst = *GeometryOptions::instance();

  std::size_t longIdx = goInst.longIdx;
  std::size_t shortIdx = goInst.shortIdx;

  std::size_t quad1Id = goInst.quad1Id;
  std::size_t quad2Id = goInst.quad2Id;
  std::size_t quad3Id = goInst.quad3Id;
  std::size_t xCorrectorId = goInst.xCorrectorId;
  std::size_t dipoleId = goInst.dipoleId;

  std::size_t nEntries = m_tree->GetEntries();
  MagneticFieldStores stores;
  stores.reserve(nEntries);
  for (std::size_t i = 0; i < nEntries; i++) {
    m_tree->GetEntry(i);

    auto store = std::make_shared<MagneticFieldStore>();
    store->store.reserve(5);

    // Quad gradients T/m
    store->store.insert(
        {quad1Id,
         Acts::MagneticFieldProvider::Cache(
             std::in_place_type<IdealQuadrupoleMagField::Cache>,
             m_quad1Grad * Acts::UnitConstants::T / Acts::UnitConstants::m)});
    store->store.insert(
        {quad2Id,
         Acts::MagneticFieldProvider::Cache(
             std::in_place_type<IdealQuadrupoleMagField::Cache>,
             m_quad2Grad * Acts::UnitConstants::T / Acts::UnitConstants::m)});
    store->store.insert(
        {quad3Id,
         Acts::MagneticFieldProvider::Cache(
             std::in_place_type<IdealQuadrupoleMagField::Cache>,
             m_quad3Grad * Acts::UnitConstants::T / Acts::UnitConstants::m)});

    // XCOR strength T
    store->store.insert(
        {dipoleId, Acts::MagneticFieldProvider::Cache(
                       std::in_place_type<ConstantMagField::Cache>,
                       m_dipoleStrength * Acts::UnitConstants::T, shortIdx)});

    // Dipole strength T
    store->store.insert(
        {xCorrectorId,
         Acts::MagneticFieldProvider::Cache(
             std::in_place_type<ConstantMagField::Cache>,
             m_xCorrectorStrength * Acts::UnitConstants::T, longIdx)});

    // Store the fields
    stores.insert({m_eventId, std::move(store)});
  }

  m_storeCollection = std::make_shared<MagneticFieldStoreCollection>(stores);
}

std::shared_ptr<MagneticFieldStoreCollection>
E320MagneticFieldParametersProvider::getMagneticFieldStoreCollection() {
  return m_storeCollection;
}

}  // namespace E320
