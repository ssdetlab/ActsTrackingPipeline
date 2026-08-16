#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

/// @brief Reads in MaterialTrack information from a root file
/// and fills it into a format to be understood by the MaterialMapping
/// algorithm
class RootMaterialTrackReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Material collection to read
    std::string outputMaterialTracks;
    /// Name of the input tree
    std::string treeName;
    /// List of input files
    std::vector<std::string> filePaths;
    /// Transform to Acts coordinate system
    Acts::Transform3 toWorldTransform;
    /// Read surface information from the root file
    bool readCachedSurfaceInformation;
  };

  /// Constructor
  ///
  /// @param config The Configuration struct
  /// @param level The log level
  RootMaterialTrackReader(const Config& config, Acts::Logging::Level level);

  /// Destructor
  ~RootMaterialTrackReader() override = default;

  /// Framework name() method
  std::string name() const override;

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const override;

  /// Read out data from the input stream
  ///
  /// @param ctx The algorithm context
  ///
  /// @return algorithm process code
  ProcessCode read(const AlgorithmContext& ctx) override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// The logger
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  /// Write data handle
  WriteDataHandle<std::vector<Acts::RecordedMaterialTrack>>
      m_outputMaterialTracks{this, "OutputMaterialTracks"};

  /// Mutex used to protect multi-threaded reads
  std::mutex m_mutex;

  /// Vector of {eventNr, entryMin, entryMax}
  std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> m_eventMap;

  /// The IO handles
  TFile* m_file = nullptr;
  TTree* m_tree = nullptr;
  TChain* m_chainOwner = nullptr;

 protected:
  /// NOTE: the naming convention here is external

  /// Event identifier.
  uint32_t m_eventId = 0;

  /// Start global x
  float m_v_x = 0;
  /// Start global y
  float m_v_y = 0;
  /// Start global z
  float m_v_z = 0;
  /// Start global momentum x
  float m_v_px = 0;
  /// Start global momentum y
  float m_v_py = 0;
  /// Start global momentum z
  float m_v_pz = 0;
  /// Start phi direction
  float m_v_phi = 0;
  /// Start eta direction
  float m_v_eta = 0;
  /// Thickness in X0
  float m_tX0 = 0;
  /// Thickness in L0
  float m_tL0 = 0;

  /// Step x position
  std::vector<float>* m_step_x = nullptr;
  /// Step y position
  std::vector<float>* m_step_y = nullptr;
  /// Step z position
  std::vector<float>* m_step_z = nullptr;
  /// Step x direction
  std::vector<float>* m_step_dx = nullptr;
  /// Step y direction
  std::vector<float>* m_step_dy = nullptr;
  /// Step z direction
  std::vector<float>* m_step_dz = nullptr;
  /// Step length
  std::vector<float>* m_step_length = nullptr;
  /// Step material X0
  std::vector<float>* m_step_X0 = nullptr;
  /// Step material L0
  std::vector<float>* m_step_L0 = nullptr;
  /// Step material A
  std::vector<float>* m_step_A = nullptr;
  /// Step material Z
  std::vector<float>* m_step_Z = nullptr;
  /// Step material rho
  std::vector<float>* m_step_rho = nullptr;

  /// ID of the surface associated with the step
  std::vector<std::uint64_t>* m_sur_id = nullptr;
  /// x position of the center of the surface associated with the step
  std::vector<float>* m_sur_x = nullptr;
  /// y position of the center of the surface associated with the step
  std::vector<float>* m_sur_y = nullptr;
  /// z position of the center of the surface associated with the step
  std::vector<float>* m_sur_z = nullptr;
  /// Path correction when associating material to the given surface
  std::vector<float>* m_sur_pathCorrection = nullptr;
};
