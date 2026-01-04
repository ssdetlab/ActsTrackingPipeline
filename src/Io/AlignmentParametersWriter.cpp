#include "TrackingPipeline/Io/AlignmentParametersWriter.hpp"

#include "Acts/Surfaces/Surface.hpp"
#include <Acts/Utilities/Logger.hpp>

#include <TVector3.h>

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

  int buf_size = 32000;
  int split_lvl = 0;

  // Detector element geometry ID
  m_tree->Branch("geometryId", &m_geoId, buf_size, split_lvl);

  // Detector element nominal transform
  m_tree->Branch("nominalTranslation", &m_nominalTranslation, buf_size,
                 split_lvl);
  m_tree->Branch("nominalRotation", &m_nominalRotation, buf_size, split_lvl);

  // Detector element new aligned transform
  m_tree->Branch("newTranslation", &m_newTranslation, buf_size, split_lvl);
  m_tree->Branch("newRotation", &m_newRotation, buf_size, split_lvl);

  // Difference between nominal and aligned
  m_tree->Branch("deltaTranslation", &m_deltaTranslation, buf_size, split_lvl);
  m_tree->Branch("deltaRotation", &m_deltaRotation, buf_size, split_lvl);

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

  ACTS_DEBUG("Received " << inputParameters.size() << " alignment results");
  std::lock_guard<std::mutex> lock(m_mutex);

  Acts::GeometryContext defGctx;
  for (const auto& [detElement, transform] : inputParameters) {
    m_geoId = detElement->surface().geometryId().sensitive();

    ACTS_DEBUG("Surface " << m_geoId);
    Acts::Vector3 nominalTrans = detElement->transform(defGctx).translation();
    Acts::RotationMatrix3 nominalRotMat =
        detElement->transform(defGctx).rotation();

    ACTS_DEBUG("Nominal translation " << nominalTrans.transpose());
    ACTS_DEBUG("Nominal rotation\n" << nominalRotMat);

    ACTS_DEBUG("Aligned translation " << transform.translation().transpose());
    ACTS_DEBUG("Aligned rotation\n" << transform.rotation());

    Acts::Vector3 deltaTrans = transform.translation() - nominalTrans;
    Acts::RotationMatrix3 deltaRotMat =
        nominalRotMat.inverse() * transform.rotation();

    ACTS_DEBUG("Delta translation " << deltaTrans.transpose());
    ACTS_DEBUG("Delta rotation\n" << deltaRotMat);

    // Store nominal
    m_nominalTranslation.SetXYZ(nominalTrans.x(), nominalTrans.y(),
                                nominalTrans.z());
    TArrayD nominalData(9);
    for (std::size_t i = 0; i < 9; i++) {
      nominalData[i] = nominalRotMat(i);
    }
    m_nominalRotation.Use(3, 3, nominalData.GetArray());

    // Store new
    m_newTranslation.SetXYZ(transform.translation().x(),
                            transform.translation().y(),
                            transform.translation().z());
    TArrayD newData(9);
    for (std::size_t i = 0; i < 9; i++) {
      newData[i] = transform.rotation()(i);
    }
    m_newRotation.Use(3, 3, newData.GetArray());

    // Store delta
    m_deltaTranslation.SetXYZ(deltaTrans.x(), deltaTrans.y(), deltaTrans.z());
    TArrayD deltaData(9);
    for (std::size_t i = 0; i < 9; i++) {
      deltaData[i] = deltaRotMat(i);
    }
    m_deltaRotation.Use(3, 3, deltaData.GetArray());

    m_tree->Fill();
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}
