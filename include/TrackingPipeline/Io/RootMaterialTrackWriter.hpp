#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <cstdint>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

class TFile;
class TTree;

namespace Acts {

// Using some short hands for Recorded Material
using RecordedMaterial = MaterialInteractor::result_type;
// And recorded material track
// - this is start:  position, start momentum
//   and the Recorded material
using RecordedMaterialTrack =
    std::pair<std::pair<Acts::Vector3, Acts::Vector3>, RecordedMaterial>;

}  // namespace Acts

/// @class RootMaterialTrackWriter
///
/// @brief Writes out MaterialTrack collections to a root file
class RootMaterialTrackWriter : public IWriter {
 public:
  struct Config {
    /// material collection to write
    std::string inputMaterialTracks = "material-tracks";
    /// path of the output file
    std::string filePath = "";
    /// file access mode
    std::string fileMode = "RECREATE";
    /// name of the output tree
    std::string treeName = "material-tracks";

    /// Re-calculate total values from individual steps (for cross-checks)
    bool recalculateTotals = false;
    /// Write aut pre and post step (for G4), otherwise central step position
    bool prePostStep = false;
    /// Write the surface to which the material step correpond
    bool storeSurface = false;
    /// Write the volume to which the material step correpond
    bool storeVolume = false;
    /// Collapse consecutive interactions of a single surface into a single
    /// interaction
    bool collapseInteractions = false;
  };

  /// Constructor with
  /// @param config configuration struct
  /// @param level logging level
  RootMaterialTrackWriter(const Config& config,
                          Acts::Logging::Level level = Acts::Logging::INFO);

  /// Virtual destructor
  ~RootMaterialTrackWriter() override;

  /// Framework initialize method
  ProcessCode finalize() override;

  std::string name() const override { return "RootMaterialTrackWriter"; };

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  ProcessCode write(const AlgorithmContext& context) override;

 private:
  /// The config class
  Config m_cfg;
  /// mutex used to protect multi-threaded writes
  std::mutex m_writeMutex;
  /// The output file name
  TFile* m_outputFile = nullptr;
  /// The output tree name
  TTree* m_outputTree = nullptr;

  /// Event identifier.
  uint32_t m_eventId = 0;

  /// start global x
  float m_v_x = 0;
  /// start global y
  float m_v_y = 0;
  /// start global z
  float m_v_z = 0;
  /// start global momentum x
  float m_v_px = 0;
  /// start global momentum y
  float m_v_py = 0;
  /// start global momentum z
  float m_v_pz = 0;
  /// start phi direction
  float m_v_phi = 0;
  /// start eta direction
  float m_v_eta = 0;
  /// thickness in X0/L0
  float m_tX0 = 0;
  /// thickness in X0/L0
  float m_tL0 = 0;

  /// step x (start) position (optional)
  std::vector<float> m_step_sx;
  /// step y (start) position (optional)
  std::vector<float> m_step_sy;
  /// step z (start) position (optional)
  std::vector<float> m_step_sz;
  /// step x position
  std::vector<float> m_step_x;
  /// step y position
  std::vector<float> m_step_y;
  /// step z position
  std::vector<float> m_step_z;
  /// step x (end) position (optional)
  std::vector<float> m_step_ex;
  /// step y (end) position (optional)
  std::vector<float> m_step_ey;
  /// step z (end) position (optional)
  std::vector<float> m_step_ez;
  /// step x direction
  std::vector<float> m_step_dx;
  /// step y direction
  std::vector<float> m_step_dy;
  /// step z direction
  std::vector<float> m_step_dz;
  /// step length
  std::vector<float> m_step_length;
  /// step material x0
  std::vector<float> m_step_X0;
  /// step material l0
  std::vector<float> m_step_L0;
  /// step material A
  std::vector<float> m_step_A;
  /// step material Z
  std::vector<float> m_step_Z;
  /// step material rho
  std::vector<float> m_step_rho;

  /// ID of the surface associated with the step
  std::vector<std::uint64_t> m_sur_id;
  /// Type of the surface associated with the step
  std::vector<int32_t> m_sur_type;
  /// x position of the center of the surface associated with the step
  std::vector<float> m_sur_x;
  /// y position of the center of the surface associated with the step
  std::vector<float> m_sur_y;
  /// z position of the center of the surface associated with the step
  std::vector<float> m_sur_z;
  /// path correction when associating material to the given surface
  std::vector<float> m_sur_pathCorrection;
  /// Min range of the surface associated with the step
  std::vector<float> m_sur_range_min;
  /// Max range of the surface associated with the step
  std::vector<float> m_sur_range_max;

  /// ID of the volume associated with the step
  std::vector<std::uint64_t> m_vol_id;

  std::unique_ptr<const Acts::Logger> m_logger;

  ReadDataHandle<std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>>
      m_inputMaterialTracks{this, "InputMaterialTracks"};
};
