#include "TrackingPipeline/Io/AlignmentParametersProvider.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>

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

  // Detector element geometry ID
  m_tree->SetBranchAddress("geometryId", &m_geoId);

  // Detector element nominal transform
  m_tree->SetBranchAddress("nominalTranslation", &m_nominalTranslation);
  m_tree->SetBranchAddress("nominalRotation", &m_nominalRotation);

  // Detector element new aligned transform
  m_tree->SetBranchAddress("newTranslation", &m_newTranslation);
  m_tree->SetBranchAddress("newRotation", &m_newRotation);

  // Number of alignment degrees of freedom
  m_tree->SetBranchAddress("alignmentDof", &m_alignmentDof);

  // Alignment parameters cov matrix
  m_tree->SetBranchAddress("alignmentCov", &m_alignmentCov);

  if (m_tree->GetEntries() != 1) {
    throw std::runtime_error(
        "Alignment tree is expected to have a single entry");
  }

  m_store = std::make_shared<AlignmentStore>();
  m_tree->GetEntry(0);

  for (std::size_t i = 0; i < m_geoId->size(); i++) {
    Acts::GeometryIdentifier geoId;
    geoId.setSensitive(m_geoId->at(i));
    if (geoId.sensitive() == 0u) {
      continue;
    }
    Acts::Vector3 translation(m_newTranslation->at(i).X(),
                              m_newTranslation->at(i).Y(),
                              m_newTranslation->at(i).Z());

    const auto& newRotation = m_newRotation->at(i);
    Acts::RotationMatrix3 rotation = Acts::RotationMatrix3::Identity();
    rotation << newRotation(0, 0), newRotation(1, 0), newRotation(2, 0),
        newRotation(0, 1), newRotation(1, 1), newRotation(2, 1),
        newRotation(0, 2), newRotation(1, 2), newRotation(2, 2);

    Acts::Transform3 transform = Acts::Transform3::Identity();
    transform.translate(translation);
    transform.rotate(rotation);

    m_store->store.emplace(geoId, transform);
  }
  m_store->covariance = Acts::ActsDynamicMatrix(m_alignmentDof, m_alignmentDof);
  for (std::size_t i = 0; i < m_alignmentDof; i++) {
    for (std::size_t j = 0; j < m_alignmentDof; j++) {
      m_store->covariance(i, j) = (*m_alignmentCov)(i, j);
    }
  }
}

std::shared_ptr<AlignmentStore>
AlignmentParametersProvider::getAlignmentStore() {
  return m_store;
}
