#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <unordered_map>

#include <RtypesCore.h>

#include "DetectorEvent.hpp"
#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TVector3.h"
#include "TrackingPipeline/EventData/ClusterMatrix.hpp"
#include "TrackingPipeline/EventData/E320EpicsParity.hpp"
#include "TrackingPipeline/EventData/E320EventMetaData.hpp"
#include "TrackingPipeline/EventData/ExtendedSourceLink.hpp"
#include "TrackingPipeline/Geometry/E320GeometryOptions.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IReader.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

namespace E320 {

/// @brief E320-specific data reader
class E320RootDataReader : public IReader {
 public:
  /// Exntended sizes shorthands
  static constexpr const std::size_t ExtendedLocalSize =
      ExtendedSourceLink::localSubspaceSize;
  static constexpr const std::size_t ExtendedGlobalSize =
      ExtendedSourceLink::globalSubspaceSize;

  /// Particle hit shorthand
  using Hit = std::pair<std::size_t, std::size_t>;

  /// @brief Cluster pixel matrixset shorthand
  using ClusterMatrixSet = std::vector<ClusterMatrix<std::uint16_t>>;

  /// @brief Measurement error model enum
  enum MeasurementErrorModel : int { Uniform = 0, Allpix2 = 1 };

  /// @brief The nested configuration struct
  struct Config {
    /// Collection with the detector measurement data
    std::string outputSourceLinks;
    /// Collection with the cluster matrices data
    std::string outputClusterMatrices;
    /// Collection with the detector measurement data indices
    std::string outputDetSourceLinkIndices;
    /// Collection with the constraint measurement data indices
    std::string outputConstraintSourceLinkIndices;
    /// Collection with the event meta data
    std::string outputEventMetaData;
    /// The names of the input files
    std::vector<std::string> filePaths;
    /// Name of the input tree
    std::string treeName;
    /// The keys we have in the ROOT file
    std::string eventKey;
    /// Epics parity requirement flag
    bool requireEpicsParity;
    /// Required Epics parity
    EpicsParity requiredEpicsParity;
    /// Measurement error model
    MeasurementErrorModel measurementErrorModel;
    /// Maximum detector occupancy
    std::size_t maxOccupancy;
    /// Geometry ID scope
    int minGeoId;
    int maxGeoId;
    /// Maximum size of a pixel cluster
    std::size_t maxClusterSize;
    /// Surface map for high-precision local to global conversion
    std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
        surfaceMap;
    /// Backwards propagation flag
    bool backwards;
  };

  E320RootDataReader(const E320RootDataReader&) = delete;
  E320RootDataReader(const E320RootDataReader&&) = delete;

  /// @brief constructor
  ///
  /// @param config configuration struct
  /// @param level logging level
  E320RootDataReader(const Config& config, Acts::Logging::Level level);

  /// @brief get reader name
  std::string name() const override { return "E320RootDataReader"; }

  /// @brief return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const final;

  /// @brief read out data from the input stream
  ///
  /// @param ctx algorithm context of the event
  ///
  /// @return execution process code
  ProcessCode read(const AlgorithmContext& ctx) override;

  /// @brief readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// @brief private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  /// WriteDataHandles for the data
  WriteDataHandle<std::vector<Acts::SourceLink>> m_outputSourceLinks{
      this, "OutputSourceLinks"};

  WriteDataHandle<ClusterMatrixSet> m_outputClusterMatrices{
      this, "OutputClusterMatrices"};

  WriteDataHandle<std::vector<std::size_t>> m_outputDetSourceLinksIndices{
      this, "OutputDetSourceLinksIndices"};

  WriteDataHandle<std::vector<std::size_t>>
      m_outputConstraintSourceLinksIndices{
          this, "OutputConstraintSourceLinksIndices"};

  WriteDataHandle<EventMetaData> m_outputEventMetaData{this, "OutputMetaData"};

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// Vector of {eventNr, entryMin, entryMax}
  std::vector<std::tuple<uint32_t, std::size_t, std::size_t>> m_eventMap;

  /// Geometry ID map between EUDAQ and ACTS
  std::unordered_map<std::uint8_t, std::size_t> m_geoIdMap;

  /// Geometry ID of the orbit reference surface
  Acts::GeometryIdentifier m_orbitReferenceSurfaceGeoId;

  /// Map of measurement cluster size to Allpix2 measurement errors
  std::unordered_map<std::size_t, std::pair<double, double>> m_allpixErrors;

  /// The input tree handle
  TTree* m_tree = nullptr;
  TFile* m_file = nullptr;
  TChain* m_chainOwner = nullptr;

  /// E320-specific geometry options
  E320::GeometryOptions m_gOpt;

 protected:
  /// Detector event
  E320Io::DetectorEvent* m_detEvent = nullptr;

  /// Estimated orbit position at the reference
  TVector3* m_orbitReferencePositionEst = nullptr;

  /// Estimated orbit momentum at the reference
  TLorentzVector* m_orbitReferenceMomentumEst = nullptr;

  /// Estimated orbit fit covariance at the reference
  TMatrixD* m_orbitReferenceCov = nullptr;

  /// Beam bunch covariance at the reference
  TMatrixD* m_orbitReferenceBeamBunchCov = nullptr;

  /// Estimated orbit covariance due to the alignment uncertainty
  TMatrixD* m_orbitReferenceAlignmentCov = nullptr;

  /// Estimated orbit covariance due to gradient-dependent mismodeling
  TMatrixD* m_orbitReferenceGradientSystematicCov = nullptr;

  /// Estimated orbit shifts for different alignment samples
  std::vector<TVectorD>* m_orbitReferenceAlignmentSampleShifts = nullptr;

  /// Estimated orbit covariance due to the alignment sample uncertainty
  TMatrixD* m_orbitReferenceAlignmentSampleCov = nullptr;

  /// Orbit fit chi2
  double m_orbitReferenceFitChi2Est = 0;

  /// Event number handle
  ULong64_t m_eventId = 0;
};

}  // namespace E320
