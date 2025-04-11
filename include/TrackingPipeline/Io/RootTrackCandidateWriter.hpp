#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IWriter.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

using namespace Acts::UnitLiterals;

/// @brief Writer to store track candidate data in
/// ROOT file
///
/// Writer that accepts track candidate data from CKF
/// derives the basic performance metrics, such as
/// chi2 and residuals, and stores them in a ROOT file.
class RootTrackCandidateWriter : public IWriter {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Surface accessor
    Acts::SourceLinkSurfaceAccessor surfaceAccessor;
    /// Track candidate collection
    std::string inputTrackCandidates;
    /// Name of the input tree
    std::string treeName;
    /// The names of the input files
    std::string filePath;
  };

  RootTrackCandidateWriter(const RootTrackCandidateWriter &) = delete;
  RootTrackCandidateWriter(const RootTrackCandidateWriter &&) = delete;

  /// @brief Constructor
  ///
  /// @param config The Configuration struct
  RootTrackCandidateWriter(const Config &config, Acts::Logging::Level level);

  /// @brief Finalize method
  ProcessCode finalize() override;

  /// Writer name() method
  std::string name() const override { return "RootTrackCandidateWriter"; }

  /// Write out data to the input stream
  ProcessCode write(const AlgorithmContext &ctx) override;

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// Private access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  ReadDataHandle<Tracks> m_inputTrackCandidates{this, "TrackCandidates"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// The output file
  TFile *m_file = nullptr;

  /// The output tree
  TTree *m_tree = nullptr;

 protected:
  /// Measurement hits
  std::vector<TVector3> m_trackHits;

  /// CKF predicted track hits
  std::vector<TVector3> m_predictedTrackHits;
  std::vector<TVector3> m_filteredTrackHits;

  /// CKF residuals with respect to the measurements
  std::vector<TVector3> m_predictedResiduals;
  std::vector<TVector3> m_filteredResiduals;

  /// CKF pulls with respect to the measurements
  std::vector<TVector3> m_predictedPulls;
  std::vector<TVector3> m_filteredPulls;

  /// Chi2 of the track
  /// with respect ot the
  /// measurement
  std::vector<double> m_chi2Predicted;
  std::vector<double> m_chi2Filtered;

  /// Number of degrees of freedom
  /// of the track
  int m_ndf;

  /// TrackId
  int m_trackId;

  /// EventId
  int m_eventId;

  /// Mutex to protect the tree filling
  std::mutex m_mutex;
};
