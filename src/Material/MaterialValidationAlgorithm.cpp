#include "TrackingPipeline/Material/MaterialValidationAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <stdexcept>

MaterialValidationAlgorithm::MaterialValidationAlgorithm(
    const MaterialValidationAlgorithm::Config& cfg, Acts::Logging::Level level)
    : IAlgorithm("MaterialValidationAlgorithm", level), m_cfg(cfg) {
  // Prepare the I/O collections
  m_inputMaterialTracks.initialize(m_cfg.inputMaterialTracks);
  m_outputMaterialTracks.initialize(m_cfg.outputMaterialTracks);

  // Check the configuration - material validater
  if (m_cfg.materialValidater == nullptr) {
    throw std::invalid_argument("Missing material validater");
  }
}

ProcessCode MaterialValidationAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const auto& mTracks = m_inputMaterialTracks(ctx);

  ACTS_DEBUG("Received " << mTracks.size() << " material tracks");

  // The output recorded material track collection
  std::vector<Acts::RecordedMaterialTrack> recordedMaterialTracks;
  recordedMaterialTracks.reserve(mTracks.size());

  // Loop over the validation tracks
  for (const auto& mTrack : mTracks) {
    const Acts::Vector3& position = mTrack.first.first;
    const Acts::Vector3& direction = mTrack.first.second;

    // Record the material
    auto rMaterial = m_cfg.materialValidater->recordMaterial(
        ctx.geoContext, ctx.magFieldContext, position, direction);
    recordedMaterialTracks.push_back(std::move(rMaterial));
  }

  ACTS_DEBUG("Sending " << recordedMaterialTracks.size() << " material tracks");

  // Send the recorded material out
  m_outputMaterialTracks(ctx, std::move(recordedMaterialTracks));

  return ProcessCode::SUCCESS;
}
