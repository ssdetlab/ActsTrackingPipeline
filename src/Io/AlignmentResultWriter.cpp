#include "TrackingPipeline/Io/AlignmentResultWriter.hpp"

#include "Acts/Surfaces/Surface.hpp"

AlignmentParametersWriter::AlignmentParametersWriter(const Config& config,
                                                     Acts::Logging::Level level)
    : m_cfg(config), m_logger(Acts::getDefaultLogger(name(), level)) {
  //  if (m_cfg.filePath.empty()) {
  //    throw std::invalid_argument("Missing filename");
  //  }
  //  if (m_cfg.treeName.empty()) {
  //    throw std::invalid_argument("Missing tree name");
  //  }

  //------------------------------------------------------------------
  // Initialize the data handles
  m_alignmentResults.initialize(m_cfg.inputAlignmentResults);
}

ProcessCode AlignmentParametersWriter::finalize() {
  if (m_file) {
    m_file->Write();
    m_file->Close();
    delete m_file;
  }
  return ProcessCode::SUCCESS;
}

ProcessCode AlignmentParametersWriter::write(const AlgorithmContext& ctx) {
  auto inputParameters = m_alignmentResults(ctx);

  std::lock_guard<std::mutex> lock(m_mutex);

  for (const auto& [detElement, transform] : inputParameters) {
    std::cout << "-------------------------------\n";
    std::cout << "Parameters: " << detElement->surface().geometryId() << "\n";
    std::cout << "Translation: " << transform.translation().transpose() << "\n";
    std::cout << "Rotation: \n" << transform.rotation() << "\n";
  }

  // Return success flag
  return ProcessCode::SUCCESS;
}

