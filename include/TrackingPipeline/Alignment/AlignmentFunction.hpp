#pragma once

/// @brief Alignment function wrapper that takes the above parameters
/// and runs the alignment procedure
#include <Acts/EventData/SourceLink.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <vector>

#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"
#include "TrackingPipeline/Utilities/NonOwningProjectionContainer.hpp"
class AlignmentFunction {
 public:
  /// Containers shorthands
  using SourceLinkContainer = std::vector<
      detail::NonOwningProjectionContainer<std::vector<Acts::SourceLink>>>;
  using TrackParametersContainer = detail::NonOwningProjectionContainer<
      std::vector<Acts::CurvilinearTrackParameters>>;
  using MagneticFieldParametersContainer = detail::NonOwningProjectionContainer<
      std::vector<std::shared_ptr<MagneticFieldStore>>>;

  /// @brief Alignment result struct
  struct AlignmentResult {
    /// Change of alignment parameters
    Acts::ActsDynamicVector deltaAlignmentParameters;
    /// Aligned parameters for detector elements
    std::unordered_map<Acts::DetectorElementBase*, Acts::Transform3>
        alignedParameters;
    /// Covariance of alignment parameters
    Acts::ActsDynamicMatrix alignmentCovariance;
    /// Average chi2/ndf (ndf is the measurement dim)
    double averageChi2ONdf = std::numeric_limits<double>::max();
    /// Total delta chi2
    double deltaChi2 = std::numeric_limits<double>::max();
    /// Final chi2
    double chi2 = 0;
    /// Measurement dimension from all tracks
    std::size_t measurementDim = 0;
    /// Alignment degree of freedom
    std::size_t alignmentDof = 0;
    /// Number of tracks used for alignment
    std::size_t numTracks = 0;
    /// Updated initial track states
    std::vector<
        Acts::GenericCurvilinearTrackParameters<Acts::ParticleHypothesis>>
        updatedInitialTrackStates;
  };

  /// @brief operator performing the alignment procedure
  ///
  /// @param gctx current geometry context
  /// @param mctx current magnetic field context
  /// @param cctx current calibration context
  /// @param trackFitSourceLinks source links participating in the alignment fit
  /// @param initialTrackStateFitSourceLinks source links participating in
  /// the initial track state fit
  /// @param initialParameters initial track parameters
  /// @param magFieldParameters track-specific magnetic field configurations
  ///
  /// @return struct containing the results of the alignment procedure
  virtual AlignmentResult operator()(
      const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
      const Acts::CalibrationContext& cctx,
      const SourceLinkContainer& alignmentFitSourceLinks,
      const SourceLinkContainer& initialTrackStateFitSourceLinks,
      const TrackParametersContainer& initialParameters,
      const MagneticFieldParametersContainer& magFieldParameters) const = 0;
};
