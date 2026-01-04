#include "TrackingPipeline/Io/AlignmentParametersProvider.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>

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

  m_store = std::make_shared<AlignmentContext::AlignmentStore>();
  for (std::size_t i = 0; i < m_tree->GetEntries(); i++) {
    m_tree->GetEntry(i);

    Acts::GeometryIdentifier geoId;
    geoId.setSensitive(m_geoId);
    if (!geoId.sensitive()) {
      continue;
    }

    Acts::Vector3 translation(m_newTranslation->X(), m_newTranslation->Y(),
                              m_newTranslation->Z());

    const auto& newRotation = *m_newRotation;
    Acts::RotationMatrix3 rotation = Acts::RotationMatrix3::Identity();
    rotation << newRotation(0, 0), newRotation(1, 0), newRotation(2, 0),
        newRotation(0, 1), newRotation(1, 1), newRotation(2, 1),
        newRotation(0, 2), newRotation(1, 2), newRotation(2, 2);

    Acts::Transform3 transform = Acts::Transform3::Identity();
    transform.translate(translation);
    transform.rotate(rotation);

    m_store->emplace(geoId, transform);
  }
}

std::shared_ptr<AlignmentContext::AlignmentStore>
AlignmentParametersProvider::getAlignmentStore() {
  return m_store;
}
