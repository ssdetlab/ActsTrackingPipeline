#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <cstdint>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

/// @brief Class writing out MaterialTrack collections to a root file
class RootMaterialTrackWriter : public IWriter {
 public:
  /// Recorded material shorthand
  using RecordedMaterial = Acts::MaterialInteractor::result_type;
  /// Material track shorthand {[originPos, originMom], RecordedMaterial}
  using RecordedMaterialTrack =
      std::pair<std::pair<Acts::Vector3, Acts::Vector3>, RecordedMaterial>;

  /// @brief nested configuration struct
  struct Config {
    /// Material collection to write
    std::string inputMaterialTracks;
    /// Output file path
    std::string filePath;
    /// Output tree name
    std::string treeName;
  };

  /// @brief Constructor
  ///
  /// @param config configuration struct
  /// @param level logging level
  RootMaterialTrackWriter(const Config& config, Acts::Logging::Level level);

  /// @brief Framework initialize method
  ProcessCode finalize() override;

  /// @brief get algorithm name
  std::string name() const override { return "RootMaterialTrackWriter"; };

  /// @brief Readonly access to the config
  const Config& config() const { return m_cfg; }

  /// @brief Framework write method
  ProcessCode write(const AlgorithmContext& ctx) override;

 private:
  /// The config class
  Config m_cfg;

  /// Mutex used to protect multi-threaded writes
  std::mutex m_mutex;

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// The output file
  TFile* m_file = nullptr;

  /// The output tree
  TTree* m_tree = nullptr;

  ReadDataHandle<std::vector<Acts::RecordedMaterialTrack>>
      m_inputMaterialTracks{this, "InputMaterialTracks"};

 protected:
  /// Event identifier
  std::size_t m_eventId = 0;

  /// Start global x
  double m_originPosX = 0;
  /// Start global y
  double m_originPosY = 0;
  /// Start global z
  double m_originPosZ = 0;
  /// Start global momentum x
  double m_originDirX = 0;
  /// Start global momentum y
  double m_originDirY = 0;
  /// Start global momentum z
  double m_originDirZ = 0;
  /// Start phi direction
  double m_originPhi = 0;
  /// Start eta direction
  double m_originTheta = 0;
  /// Total track thickness in X0
  double m_thicknessX0 = 0;

  /// Step x position
  std::vector<double> m_stepPosX;
  /// Step y position
  std::vector<double> m_stepPosY;
  /// Step z position
  std::vector<double> m_stepPosZ;
  /// Step x direction
  std::vector<double> m_stepDirX;
  /// Step y direction
  std::vector<double> m_stepDirY;
  /// Step z direction
  std::vector<double> m_stepDirZ;
  /// Step length
  std::vector<double> m_stepLength;
  /// Step material X0
  std::vector<double> m_stepX0;
  /// Step thickness X0
  std::vector<double> m_stepThicknessX0;
  /// Step material A
  std::vector<double> m_stepA;
  /// Step material Z
  std::vector<double> m_stepZ;
  /// Step material rho
  std::vector<double> m_stepRho;

  /// ID of the surface associated with the step
  std::vector<std::uint64_t> m_surfaceId;
  /// Type of the surface associated with the step
  std::vector<int32_t> m_surfaceType;
  /// x position of the center of the surface associated with the step
  std::vector<double> m_surfaceX;
  /// y position of the center of the surface associated with the step
  std::vector<double> m_surfaceY;
  /// z position of the center of the surface associated with the step
  std::vector<double> m_surfaceZ;
  /// Path correction when associating material to the given surface
  std::vector<double> m_surfacePathCorrection;
};
