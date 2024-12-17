#include "TrackingPipeline/Material/CoreMaterialMapping.hpp"
#include "TrackingPipeline/Material/IMaterialWriter.hpp"

#include <stdexcept>
#include <unordered_map>

CoreMaterialMapping::CoreMaterialMapping(
    const CoreMaterialMapping::Config& cfg,
    Acts::Logging::Level level)
        : IAlgorithm("CoreMaterialMapping", level), m_cfg(cfg) {
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

CoreMaterialMapping::~CoreMaterialMapping() {
    Acts::DetectorMaterialMaps detectorMaterial =
        m_cfg.materialMapper->finalizeMaps(*m_mappingState);
    // Loop over the available writers and write the maps
    for (auto& imw : m_cfg.materiaMaplWriters) {
        imw->writeMaterial(detectorMaterial);
    }
}

ProcessCode CoreMaterialMapping::execute(
    const AlgorithmContext& context) const {
        // Take the collection from the EventStore: input collection
        std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>
            mtrackCollection = m_inputMaterialTracks(context);
        
        // Write the output collections to the Event store : mapped and unmapped
        std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>
            mappedTrackCollection;
        
        std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>
            unmappedTrackCollection;
        
        // To make it work with the framework needs a lock guard
        auto mappingState =
            const_cast<Acts::MaterialMapper::State*>(m_mappingState.get());
        
        for (auto& [idTrack, mTrack] : mtrackCollection) {
            auto [mapped, unmapped] = m_cfg.materialMapper->mapMaterial(
                *mappingState, context.geoContext, context.magFieldContext, mTrack);

            mappedTrackCollection.emplace_hint(
                mappedTrackCollection.end(), 
                idTrack,
                mapped);
            unmappedTrackCollection.emplace_hint(
                unmappedTrackCollection.end(), 
                idTrack,
                unmapped);
        }
        
        // Write the mapped and unmapped material tracks to the output
        m_outputMappedMaterialTracks(context, std::move(mappedTrackCollection));
        m_outputUnmappedMaterialTracks(context, std::move(unmappedTrackCollection));
        
        return ProcessCode::SUCCESS;
}
