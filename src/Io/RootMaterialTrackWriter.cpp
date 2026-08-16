#include "TrackingPipeline/Io/RootMaterialTrackWriter.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <mutex>
#include <stdexcept>

RootMaterialTrackWriter::RootMaterialTrackWriter(
    const RootMaterialTrackWriter::Config& config, Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("RootMaterialTrackWriter", level)) {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_file = new TFile(m_cfg.filePath.c_str(), "RECREATE");
  m_tree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());

  //------------------------------------------------------------------
  // Set tree branches

  int bufSize = 32000;
  int splitLvl = 0;

  // Event identifier
  m_tree->Branch("eventId", &m_eventId, bufSize, splitLvl);

  // Start global x
  m_tree->Branch("originPosX", &m_originPosX, bufSize, splitLvl);
  // Start global y
  m_tree->Branch("originPosY", &m_originPosY, bufSize, splitLvl);
  // Start global z
  m_tree->Branch("originPosZ", &m_originPosZ, bufSize, splitLvl);
  // Start global momentum x
  m_tree->Branch("originDirX", &m_originDirX, bufSize, splitLvl);
  // Start global momentum y
  m_tree->Branch("originDirY", &m_originDirY, bufSize, splitLvl);
  // Start global momentum z
  m_tree->Branch("originDirZ", &m_originDirZ, bufSize, splitLvl);
  // Start phi direction
  m_tree->Branch("originPhi", &m_originPhi, bufSize, splitLvl);
  // Start eta direction
  m_tree->Branch("originTheta", &m_originTheta, bufSize, splitLvl);
  // Thickness in X0
  m_tree->Branch("thicknessX0", &m_thicknessX0, bufSize, splitLvl);

  // Step x position
  m_tree->Branch("stepPosX", &m_stepPosX, bufSize, splitLvl);
  // Step y position
  m_tree->Branch("stepPosY", &m_stepPosY, bufSize, splitLvl);
  // Step z position
  m_tree->Branch("stepPosZ", &m_stepPosZ, bufSize, splitLvl);
  /// Step x direction
  m_tree->Branch("stepDirX", &m_stepDirX, bufSize, splitLvl);
  // Step y direction
  m_tree->Branch("stepDirY", &m_stepDirY, bufSize, splitLvl);
  // Step z direction
  m_tree->Branch("stepDirZ", &m_stepDirZ, bufSize, splitLvl);
  // Step length
  m_tree->Branch("stepLength", &m_stepLength, bufSize, splitLvl);
  // Step material X0
  m_tree->Branch("stepX0", &m_stepX0, bufSize, splitLvl);
  // Step thickness X0
  m_tree->Branch("stepThicknessX0", &m_stepThicknessX0, bufSize, splitLvl);
  // Step material A
  m_tree->Branch("stepA", &m_stepA, bufSize, splitLvl);
  // Step material Z
  m_tree->Branch("stepZ", &m_stepZ, bufSize, splitLvl);
  // Step material rho
  m_tree->Branch("stepRho", &m_stepRho, bufSize, splitLvl);

  // ID of the surface associated with the step
  m_tree->Branch("surfaceId", &m_surfaceId, bufSize, splitLvl);
  // Type of the surface associated with the step
  m_tree->Branch("surfaceType", &m_surfaceType, bufSize, splitLvl);
  // x position of the center of the surface associated with the step
  m_tree->Branch("surfaceX", &m_surfaceX, bufSize, splitLvl);
  // y position of the center of the surface associated with the step
  m_tree->Branch("surfaceY", &m_surfaceY, bufSize, splitLvl);
  // z position of the center of the surface associated with the step
  m_tree->Branch("surfaceZ", &m_surfaceZ, bufSize, splitLvl);
  // Path correction when associating material to the given surface
  m_tree->Branch("surfacePathCorrection", &m_surfacePathCorrection, bufSize,
                 splitLvl);

  //------------------------------------------------------------------
  // Initialize the data handles
  m_inputMaterialTracks.initialize(m_cfg.inputMaterialTracks);
}

ProcessCode RootMaterialTrackWriter::finalize() {
  if (m_file != nullptr) {
    m_file->Write();
    m_file->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootMaterialTrackWriter::write(const AlgorithmContext& ctx) {
  const auto& materialTracks = m_inputMaterialTracks(ctx);

  std::lock_guard<std::mutex> lock(m_mutex);

  m_eventId = ctx.eventNumber;

  // Loop over the material tracks and write them out
  for (const auto& [originInfo, recordedMaterial] : materialTracks) {
    // Step x position
    m_stepPosX.clear();
    // Step y position
    m_stepPosY.clear();
    // Step z position
    m_stepPosZ.clear();
    // Step x direction
    m_stepDirX.clear();
    // Step y direction
    m_stepDirY.clear();
    // Step z direction
    m_stepDirZ.clear();
    // Step length
    m_stepLength.clear();
    // Step material X0
    m_stepX0.clear();
    // Step thickness X0
    m_stepThicknessX0.clear();
    // Step material A
    m_stepA.clear();
    // Step material Z
    m_stepZ.clear();
    // Step material rho
    m_stepRho.clear();

    // ID of the surface associated with the step
    m_surfaceId.clear();
    // Type of the surface associated with the step
    m_surfaceType.clear();
    // x position of the center of the surface associated with the step
    m_surfaceX.clear();
    // y position of the center of the surface associated with the step
    m_surfaceY.clear();
    // z position of the center of the surface associated with the step
    m_surfaceZ.clear();
    // Path correction when associating material to the given surface
    m_surfacePathCorrection.clear();

    // Reserve the vectors
    const auto& materialInteractions = recordedMaterial.materialInteractions;
    std::size_t nInteractions = materialInteractions.size();

    // Step x position
    m_stepPosX.reserve(nInteractions);
    // Step y position
    m_stepPosY.reserve(nInteractions);
    // Step z position
    m_stepPosZ.reserve(nInteractions);
    // Step x direction
    m_stepDirX.reserve(nInteractions);
    // Step y direction
    m_stepDirY.reserve(nInteractions);
    // Step z direction
    m_stepDirZ.reserve(nInteractions);
    // Step length
    m_stepLength.reserve(nInteractions);
    // Step material X0
    m_stepX0.reserve(nInteractions);
    // Step thickness X0
    m_stepThicknessX0.reserve(nInteractions);
    // Step material A
    m_stepA.reserve(nInteractions);
    // Step material Z
    m_stepZ.reserve(nInteractions);
    // Step material rho
    m_stepRho.reserve(nInteractions);

    // ID of the surface associated with the step
    m_surfaceId.reserve(nInteractions);
    // Type of the surface associated with the step
    m_surfaceType.reserve(nInteractions);
    // x position of the center of the surface associated with the step
    m_surfaceX.reserve(nInteractions);
    // y position of the center of the surface associated with the step
    m_surfaceY.reserve(nInteractions);
    // z position of the center of the surface associated with the step
    m_surfaceZ.reserve(nInteractions);
    // Path correction when associating material to the given surface
    m_surfacePathCorrection.reserve(nInteractions);

    // Total track thickness in X0
    m_thicknessX0 = recordedMaterial.materialInX0;

    // Set the track information at vertex
    const auto& [originPos, originDir] = originInfo;
    m_originPosX = originPos.x();
    m_originPosY = originPos.y();
    m_originPosZ = originPos.z();
    m_originDirX = originDir.x();
    m_originDirY = originDir.y();
    m_originDirZ = originDir.z();
    m_originPhi = Acts::VectorHelpers::phi(originDir);
    m_originTheta = Acts::VectorHelpers::theta(originDir);

    // Loop over the material
    for (const auto& mInt : materialInteractions) {
      const Acts::Vector3& position = mInt.position;
      const Acts::Vector3& direction = mInt.direction;

      // The material step position information
      m_stepPosX.push_back(position.x());
      m_stepPosY.push_back(position.y());
      m_stepPosZ.push_back(position.z());
      m_stepDirX.push_back(direction.x());
      m_stepDirY.push_back(direction.y());
      m_stepDirZ.push_back(direction.z());

      // Store surface information
      const Acts::Surface* surface = mInt.surface;
      if (mInt.intersectionID.value() != 0) {
        if (mInt.intersectionID.sensitive() != 0u) {
          m_surfaceId.push_back(mInt.intersectionID.sensitive());
        } else {
          m_surfaceId.push_back(mInt.intersectionID.passive());
        }
        m_surfacePathCorrection.push_back(mInt.pathCorrection);
        m_surfaceX.push_back(mInt.intersection.x());
        m_surfaceY.push_back(mInt.intersection.y());
        m_surfaceZ.push_back(mInt.intersection.z());
      } else if (surface != nullptr) {
        auto sfIntersection =
            surface
                ->intersect(
                    ctx.geoContext, mInt.position, mInt.direction,
                    Acts::BoundaryTolerance::AbsoluteCartesian(1e-4, 1e-4))
                .closest();
        m_surfaceId.push_back(surface->geometryId().value());
        m_surfacePathCorrection.push_back(1.0);
        m_surfaceX.push_back(sfIntersection.position().x());
        m_surfaceY.push_back(sfIntersection.position().y());
        m_surfaceZ.push_back(sfIntersection.position().z());
      } else {
        m_surfaceId.push_back(Acts::GeometryIdentifier().value());
        m_surfaceX.push_back(0);
        m_surfaceY.push_back(0);
        m_surfaceZ.push_back(0);
        m_surfacePathCorrection.push_back(1.0);
      }
      if (surface != nullptr) {
        m_surfaceType.push_back(surface->type());
      } else {
        m_surfaceType.push_back(-1);
      }

      // The material information
      const auto& mProps = mInt.materialSlab;
      m_stepLength.push_back(mProps.thickness());
      m_stepX0.push_back(mProps.material().X0());
      m_stepThicknessX0.push_back(mProps.thicknessInX0());
      m_stepA.push_back(mProps.material().Ar());
      m_stepZ.push_back(mProps.material().Z());
      m_stepRho.push_back(mProps.material().massDensity());
    }
    m_tree->Fill();
  }

  // return success
  return ProcessCode::SUCCESS;
}
