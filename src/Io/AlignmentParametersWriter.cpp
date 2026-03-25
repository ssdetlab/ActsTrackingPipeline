#include "TrackingPipeline/Io/AlignmentParametersWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>

#include "TVector3.h"

AlignmentParametersWriter::AlignmentParametersWriter(const Config& config,
                                                     Acts::Logging::Level level)
    : m_cfg(config), m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_file = new TFile(m_cfg.filePath.c_str(), "RECREATE");
  m_tree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());

  int bufSize = 32000;
  int splitLvl = 0;

  // Detector element geometry ID
  m_tree->Branch("geometryId", &m_geoId, bufSize, splitLvl);

  // Detector element nominal transform
  m_tree->Branch("nominalTranslation", &m_nominalTranslation, bufSize,
                 splitLvl);
  m_tree->Branch("nominalRotation", &m_nominalRotation, bufSize, splitLvl);

  // Detector element new aligned transform
  m_tree->Branch("newTranslation", &m_newTranslation, bufSize, splitLvl);
  m_tree->Branch("newRotation", &m_newRotation, bufSize, splitLvl);

  // Number of alignment degrees of freedom
  m_tree->Branch("alignmentDof", &m_alignmentDof, bufSize, splitLvl);

  // Alignment parameters cov matrix
  m_tree->Branch("alignmentCov", &m_alignmentCov, bufSize, splitLvl);

  //------------------------------------------------------------------
  // Initialize the data handles
  m_alignmentResults.initialize(m_cfg.inputAlignmentResults);
}

ProcessCode AlignmentParametersWriter::finalize() {
  m_file->Write();
  m_file->Close();
  return ProcessCode::SUCCESS;
}

ProcessCode AlignmentParametersWriter::write(const AlgorithmContext& ctx) {
  auto inputParameters = m_alignmentResults(ctx);

  ACTS_DEBUG("Received " << inputParameters.alignedParameters.size()
                         << " alignment results");
  std::lock_guard<std::mutex> lock(m_mutex);

  Acts::GeometryContext defGctx;
  for (const auto& [detElement, transform] :
       inputParameters.alignedParameters) {
    m_geoId.push_back(detElement->surface().geometryId().sensitive());

    Acts::Vector3 nominalTranslation =
        detElement->transform(defGctx).translation();
    Acts::RotationMatrix3 nominalRotation =
        detElement->transform(defGctx).rotation();

    // Store nominal
    m_nominalTranslation.push_back(TVector3(nominalTranslation.x(),
                                            nominalTranslation.y(),
                                            nominalTranslation.z()));

    TMatrixD nominalRot;
    TArrayD nominalRotationData(9);
    for (std::size_t i = 0; i < 9; i++) {
      nominalRotationData[i] = nominalRotation(i);
    }
    nominalRot.Use(3, 3, nominalRotationData.GetArray());
    m_nominalRotation.push_back(nominalRot);

    // Store new
    Acts::Vector3 newTranslation = transform.translation();
    Acts::RotationMatrix3 newRotation = transform.rotation();
    m_newTranslation.push_back(
        TVector3(newTranslation.x(), newTranslation.y(), newTranslation.z()));

    TMatrixD newRot;
    TArrayD newData(9);
    for (std::size_t i = 0; i < 9; i++) {
      newData[i] = newRotation(i);
    }
    newRot.Use(3, 3, newData.GetArray());
    m_newRotation.push_back(newRot);
  }

  m_alignmentDof = inputParameters.alignmentDof;
  const auto& alignmentCov = inputParameters.alignmentCovariance;
  TArrayD alignmentCovData(m_alignmentDof * m_alignmentDof);
  for (std::size_t i = 0; i < m_alignmentDof * m_alignmentDof; i++) {
    alignmentCovData[i] = alignmentCov(i);
  }
  m_alignmentCov.Use(m_alignmentDof, m_alignmentDof,
                     alignmentCovData.GetArray());

  m_tree->Fill();

  // Return success flag
  return ProcessCode::SUCCESS;
}
