#include "TrackingPipeline/Material/MaterialValidation.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <stdexcept>

MaterialValidation::MaterialValidation(const MaterialValidation::Config& cfg,
                                       Acts::Logging::Level level)
    : IAlgorithm("MaterialValidation", level), m_cfg(cfg) {
  // Prepare the I/O collections
  m_inputMaterialTracks.initialize(m_cfg.inputMaterialTracks);
  m_outputMaterialTracks.initialize(m_cfg.outputMaterialTracks);
  // Check the configuration - material validater
  if (m_cfg.materialValidater == nullptr) {
    throw std::invalid_argument("Missing material validater.");
  }
}

ProcessCode MaterialValidation::execute(const AlgorithmContext& context) const {
  // Take the collection from the EventStore: input collection
  std::unordered_map<std::size_t, Acts::RecordedMaterialTrack> mtracks =
      m_inputMaterialTracks(context);

  // The output recorded material track collection
  std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>
      recordedMaterialTracks;

  // Loop over the validation tracks
  for (auto& [idTrack, mTrack] : mtracks) {
    Acts::Vector3 position = mTrack.first.first;
    Acts::Vector3 direction = mTrack.first.second;

    // Record the material
    auto rMaterial = m_cfg.materialValidater->recordMaterial(
        context.geoContext, context.magFieldContext, position, direction);

    recordedMaterialTracks.emplace_hint(recordedMaterialTracks.end(), idTrack,
                                        rMaterial);
  }

  // Write the mapped and unmapped material tracks to the output
  m_outputMaterialTracks(context, std::move(recordedMaterialTracks));

  return ProcessCode::SUCCESS;
}
