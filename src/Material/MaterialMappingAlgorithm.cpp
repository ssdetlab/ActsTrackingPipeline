#include "TrackingPipeline/Material/MaterialMappingAlgorithm.hpp"

#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <stdexcept>

#include "TrackingPipeline/Material/IMaterialWriter.hpp"

MaterialMappingAlgorithm::MaterialMappingAlgorithm(
    const MaterialMappingAlgorithm::Config& cfg, Acts::Logging::Level level)
    : IAlgorithm("MaterialMappingAlgorithm", level), m_cfg(cfg) {
  // Prepare the I/O collections
  m_inputMaterialTracks.initialize(m_cfg.inputMaterialTracks);
  m_outputMappedMaterialTracks.initialize(m_cfg.mappedMaterialTracks);
  m_outputUnmappedMaterialTracks.initialize(m_cfg.unmappedMaterialTracks);

  ACTS_INFO("This algorithm requires inter-event information, "
            << "run in single-threaded mode!");

  if (m_cfg.materialMapper == nullptr) {
    throw std::invalid_argument("Missing material mapper");
  }
  // Create the state object
  m_mappingState = m_cfg.materialMapper->createState();
}

MaterialMappingAlgorithm::~MaterialMappingAlgorithm() {
  Acts::DetectorMaterialMaps detectorMaterial =
      m_cfg.materialMapper->finalizeMaps(*m_mappingState);
  // Loop over the available writers and write the maps
  for (auto& imw : m_cfg.materiaMaplWriters) {
    imw->writeMaterial(detectorMaterial);
  }
}

ProcessCode MaterialMappingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const auto& mTrackCollection = m_inputMaterialTracks(ctx);

  ACTS_DEBUG("Received " << mTrackCollection.size() << " material tracks");

  // Output collections
  std::vector<Acts::RecordedMaterialTrack> mappedTrackCollection;
  mappedTrackCollection.reserve(mTrackCollection.size());

  std::vector<Acts::RecordedMaterialTrack> unmappedTrackCollection;
  unmappedTrackCollection.reserve(mTrackCollection.size());

  // To make it work with the framework needs a lock guard
  auto* mappingState =
      const_cast<Acts::MaterialMapper::State*>(m_mappingState.get());

  std::size_t mappedInteractions = 0;
  std::size_t unmappedInteractions = 0;
  for (const auto& mTrack : mTrackCollection) {
    auto [mapped, unmapped] = m_cfg.materialMapper->mapMaterial(
        *mappingState, ctx.geoContext, ctx.magFieldContext, mTrack,
        m_cfg.materialMapperOptions);

    mappedInteractions += mapped.second.materialInteractions.size();
    unmappedInteractions += unmapped.second.materialInteractions.size();

    mappedTrackCollection.push_back(std::move(mapped));
    unmappedTrackCollection.push_back(std::move(unmapped));
  }

  ACTS_DEBUG("Sending " << mappedInteractions
                        << " mapped material interactions");
  ACTS_DEBUG("Sending " << unmappedInteractions
                        << " unmapped material interactions");

  // Write the mapped and unmapped material tracks to the output
  m_outputMappedMaterialTracks(ctx, std::move(mappedTrackCollection));
  m_outputUnmappedMaterialTracks(ctx, std::move(unmappedTrackCollection));

  return ProcessCode::SUCCESS;
}
