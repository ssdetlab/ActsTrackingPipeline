#include "TrackingPipeline/Io/AlignmentParametersProvider.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>

#include "TVector3.h"

AlignmentParametersProvider::AlignmentParametersProvider(const Config& config)
    : m_cfg(config) {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_file = new TFile(m_cfg.filePath.c_str(), "READ");
  m_tree = m_file->Get<TTree>(m_cfg.treeName.c_str());

  int buf_size = 32000;
  int split_lvl = 0;

  // Detector element geometry ID
  m_tree->SetBranchAddress("geometryId", &m_geoId);

  // Detector element nominal transform
  m_tree->SetBranchAddress("nominalTranslation", &m_nominalTranslation);
  m_tree->SetBranchAddress("nominalRotation", &m_nominalRotation);

  // Detector element new aligned transform
  m_tree->SetBranchAddress("newTranslation", &m_newTranslation);
  m_tree->SetBranchAddress("newRotation", &m_newRotation);

  // Difference between nominal and aligned
  m_tree->SetBranchAddress("deltaTranslation", &m_deltaTranslation);
  m_tree->SetBranchAddress("deltaRotation", &m_deltaRotation);

  m_store = std::make_shared<AlignmentParameters>();
  for (std::size_t i = 0; i < m_tree->GetEntries(); i++) {
    m_tree->GetEntry(i);

    Acts::GeometryIdentifier geoId;
    geoId.setSensitive(m_geoId);

    Acts::Vector3 deltaTranslation(m_deltaTranslation->X(),
                                   m_deltaTranslation->Y(),
                                   m_deltaTranslation->Z());

    Acts::RotationMatrix3 deltaRotation = Acts::RotationMatrix3::Identity();
    deltaRotation << (*m_deltaRotation)(0, 0), (*m_deltaRotation)(1, 0),
        (*m_deltaRotation)(2, 0), (*m_deltaRotation)(0, 1),
        (*m_deltaRotation)(1, 1), (*m_deltaRotation)(2, 1),
        (*m_deltaRotation)(0, 2), (*m_deltaRotation)(1, 2),
        (*m_deltaRotation)(2, 2);

    m_store->emplace(geoId, std::make_pair(deltaTranslation, deltaRotation));
  }
}

std::pair<Acts::Vector3, Acts::RotationMatrix3>&
AlignmentParametersProvider::getAlignedTransform(
    const Acts::GeometryIdentifier& geoId) {
  if (!m_store->contains(geoId)) {
    throw std::runtime_error("No detector element found in store");
  }
  return m_store->at(geoId);
}
